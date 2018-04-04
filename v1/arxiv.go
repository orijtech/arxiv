// Copyright 2017 orijtech. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package arxiv

import (
	"context"
	"encoding/xml"
	"errors"
	"fmt"
	"io/ioutil"
	"net/http"
	"sync"
	"time"

	"golang.org/x/tools/blog/atom"

	"go.opencensus.io/plugin/ochttp"
	"go.opencensus.io/trace"

	"github.com/orijtech/otils"
)

type Client struct {
	sync.RWMutex

	rt http.RoundTripper
}

func (c *Client) httpClient() *http.Client {
	c.RLock()
	defer c.RUnlock()

	rt := c.rt
	if rt == nil {
		rt = http.DefaultTransport
	}
	return &http.Client{Transport: &ochttp.Transport{Base: rt}}
}

func (c *Client) SetHTTPRoundTripper(rt http.RoundTripper) {
	c.Lock()
	defer c.Unlock()

	c.rt = rt
}

const (
	baseURL = "http://export.arxiv.org/api"
)

const recommendedThrottleDuration = time.Duration(3 * time.Second)

type SortOrder string

const (
	SortAscending  SortOrder = "ascending"
	SortDescending SortOrder = "descending"
)

type SortBy string

const (
	SortByRelevance       SortBy = "relevance"
	SortBySubmittedDate   SortBy = "submittedDate"
	SortByLastUpdatedDate SortBy = "lastUpdatedDate"
)

type Query struct {
	MaxPageNumber int64 `json:"max_page"`

	Terms string `json:"search_query"`

	Filters Filters `json:"filters"`

	PageNumber int64 `json:"start"`

	MaxResultsPerPage int64 `json:"max_results"`

	// ArticleIDs is the specific set of articles
	// to be searched for by the API.
	ArticleIDs []string `json:"id_list"`

	SortBy    SortBy    `json:"sortBy"`
	SortOrder SortOrder `json:"sortOrder"`

	// ThrottleSeconds determines the period that it'll
	// wait before making the next pagination request.
	// arxiv.org recommends that you set it to 3 seconds
	// to play nice. If this value is unset, the recommended
	// 3 second throttle duration will be used.
	// Set it to a negative value e.g -1 for no throttling.
	ThrottleSeconds int64 `json:"throttle_seconds"`
}

// Fresh struct here to avoid sending
// unnecessary query string values from Query
type reqPage struct {
	Terms      string   `json:"search_query"`
	Start      int64    `json:"start"`
	MaxResults int64    `json:"max_results"`
	ArticleIDs []string `json:"id_list"`

	SortBy    SortBy    `json:"sortBy"`
	SortOrder SortOrder `json:"sortOrder"`
}

func (q *Query) reqPage() *reqPage {
	return &reqPage{
		Terms:      q.Terms,
		Start:      q.PageNumber,
		SortBy:     q.SortBy,
		MaxResults: q.MaxResultsPerPage,
		ArticleIDs: q.ArticleIDs[:],
		SortOrder:  q.SortOrder,
	}
}

var (
	errNilQuery            = errors.New("non-nil query expected")
	errExpectingIDsOrQuery = errors.New("expecting either ArticleIDs or Terms or Filters to have been set")
)

func (q *Query) Validate() error {
	if q == nil {
		return errNilQuery
	}
	if len(q.ArticleIDs) == 0 && q.Terms == "" && len(q.Filters) == 0 {
		return errExpectingIDsOrQuery
	}
	return nil
}

type ResultsPage struct {
	Feed *atom.Feed `json:"feed"`

	PageNumber int64 `json:"page_number"`

	Err error `json:"error"`
}

const (
	defaultMaxResultsPerPage = 10
)

func (q *Query) prepareForPagination() {
	// Pages are zero-based so anything
	// less than that should be set to zero
	if q.PageNumber <= 0 {
		q.PageNumber = 0
	}

	if q.MaxResultsPerPage <= 0 {
		q.MaxResultsPerPage = defaultMaxResultsPerPage
	}
}

type cancelFn func()

var defaultClient = new(Client)

func Search(ctx context.Context, q *Query) (chan *ResultsPage, cancelFn, error) {
	return defaultClient.Search(ctx, q)
}

func (c *Client) Search(ctx context.Context, q *Query) (chan *ResultsPage, cancelFn, error) {
	ctx, span := trace.StartSpan(ctx, "arxiv/v1.Search")
	defer span.End()

	if err := q.Validate(); err != nil {
		return nil, nil, err
	}

	maxPageNumber := q.MaxPageNumber
	maxPageExceeded := func(pageNumber int64) bool {
		return maxPageNumber > 0 && pageNumber >= maxPageNumber
	}

	responsesChan := make(chan *ResultsPage)
	cancelFn, cancelChan := makeCanceler()

	var throttleDuration = recommendedThrottleDuration
	secs := q.ThrottleSeconds
	if secs >= 1 {
		throttleDuration = time.Duration(secs) * time.Second
	} else if secs <= -1 {
		throttleDuration = 0 * time.Second
	}

	// Take out the filters on their own
	filters := q.Filters[:]
	filtersStr := filters.String()

	go func() {
		defer close(responsesChan)

		pager := q
		pager.prepareForPagination()

		for {
			respPage := new(ResultsPage)
			respPage.PageNumber = pager.PageNumber
			qv, err := otils.ToURLValues(pager.reqPage())
			if err != nil {
				respPage.Err = err
				responsesChan <- respPage
				return
			}
			if len(filters) > 0 {
				qv.Set("search_query", filtersStr)
			}

			cctx, span := trace.StartSpan(ctx, "paging")
			span.Annotate([]trace.Attribute{
				trace.Int64Attribute("page", pager.PageNumber),
			}, "paging")
			fullURL := fmt.Sprintf("%s/query?%s", baseURL, qv.Encode())
			req, err := http.NewRequest("GET", fullURL, nil)
			if err != nil {
				span.End()
				respPage.Err = err
				responsesChan <- respPage
				return
			}

			slurp, _, err := c.doReqAndSlurp(cctx, req)
			span.End()
			if err != nil {
				respPage.Err = err
				responsesChan <- respPage
				return
			}

			_, umSpan := trace.StartSpan(ctx, "xml_unmarshal-to-atom.Feed")
			feed := new(atom.Feed)
			err = xml.Unmarshal(slurp, feed)
			umSpan.End()

			if err != nil {
				respPage.Err = err
				responsesChan <- respPage
				return
			}

			if feed == nil || len(feed.Entry) == 0 {
				// No more content since pages must be contiguous
				// before we encounter the first empty page.
				break
			}

			respPage.Feed = feed
			responsesChan <- respPage

			pager.PageNumber += 1
			if maxPageExceeded(pager.PageNumber) {
				break
			}

			select {
			case <-time.After(throttleDuration):
			case <-cancelChan:
				return
			}
		}
	}()

	return responsesChan, cancelFn, nil
}

func makeCanceler() (cancelFn, <-chan bool) {
	cancelChan := make(chan bool)
	var cancelOnce sync.Once
	cancelFn := func() { cancelOnce.Do(func() { close(cancelChan) }) }
	return cancelFn, cancelChan
}

func (c *Client) doReqAndSlurp(ctx context.Context, req *http.Request) ([]byte, http.Header, error) {
	ctx, span := trace.StartSpan(ctx, "arxiv/v1/doReqAndslurp")
	defer span.End()

	res, err := c.httpClient().Do(req)
	if err != nil {
		return nil, nil, err
	}

	if res.Body != nil {
		defer res.Body.Close()
	}

	if !otils.StatusOK(res.StatusCode) {
		msg := res.Status
		if res.Body != nil {
			slurp, _ := ioutil.ReadAll(res.Body)
			if len(slurp) > 4 {
				msg = string(slurp)
			}
		}
		return nil, res.Header, errors.New(msg)
	}

	slurp, err := ioutil.ReadAll(res.Body)
	return slurp, res.Header, err
}
