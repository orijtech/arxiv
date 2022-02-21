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

package arxiv_test

import (
	"bytes"
	"encoding/xml"
	"fmt"
	"io/ioutil"
	"net/http"
	"os"
	"strconv"
	"testing"
)

func TestSearch(t *testing.T) {
	client := new(arxiv.Client)
	backend := &testBackend{route: searchRoute}
	client.SetHTTPRoundTripper(backend)

	tests := [...]struct {
		req       *arxiv.Query
		wantErr   bool
		wantPages []*arxiv.ResultsPage
	}{
		0: {
			req: &arxiv.Query{
				Terms:         "deep dream",
				MaxPageNumber: 3,

				ThrottleSeconds: -1,
			},
			wantPages: []*arxiv.ResultsPage{
				resultsPageFromFile("deep dream", 0),
				resultsPageFromFile("deep dream", 1),
				resultsPageFromFile("deep dream", 2),
			},
		},
		1: {
			req:     nil,
			wantErr: true,
		},
		2: {
			req:     &arxiv.Query{},
			wantErr: true,
		},
	}

	for i, tt := range tests {
		resPageChan, cancel, err := client.Search(tt.req)
		if tt.wantErr {
			if err == nil {
				t.Errorf("#%d want non-nil error", i)
			}
			continue
		}

		if err != nil {
			t.Errorf("#%d:: got err: %v", i, err)
			continue
		}

		if cancel == nil {
			t.Errorf("#%d:: expecting a non-nil canceler", i)
			continue
		}

		var gotPages []*arxiv.ResultsPage
		for resPage := range resPageChan {
			gotPages = append(gotPages, resPage)
		}
		cancel()

		gotBlob, wantBlob := blobify(gotPages), blobify(tt.wantPages)
		if !bytes.Equal(gotBlob, wantBlob) {
			t.Errorf("#%d:\ngotBlob:  %s\nwantBlob: %s\n", i, gotBlob, wantBlob)
		}
	}
}

func blobify(v interface{}) []byte {
	blob, _ := xml.Marshal(v)
	return blob
}

const (
	searchRoute = "search"
)

type testBackend struct {
	route string
}

var _ http.RoundTripper = (*testBackend)(nil)

func (tb *testBackend) RoundTrip(req *http.Request) (*http.Response, error) {
	switch tb.route {
	case searchRoute:
		return tb.searchRoundTrip(req)
	default:
		return makeResp("unknown route", http.StatusNotFound), nil
	}
}

func makeResp(status string, code int) *http.Response {
	return &http.Response{
		StatusCode: code,
		Status:     status,
		Body:       nil,
	}
}

func (tb *testBackend) searchRoundTrip(req *http.Request) (*http.Response, error) {
	query := req.URL.Query()
	terms := query.Get("search_query")
	if terms == "" {
		return makeResp(`expecting "search_query"`, http.StatusBadRequest), nil
	}

	pageNumber, err := strconv.ParseInt(query.Get("start"), 10, 64)
	if err != nil {
		msg := fmt.Sprintf(`trying to parse "start" err: %q. Expecting an integer`, err)
		return makeResp(msg, http.StatusBadRequest), nil
	}

	diskPath := resultsPagePath(terms, pageNumber)
	return makeRespFromFile(diskPath), nil
}

func resultsPagePath(searchTerm string, pageNumber int64) string {
	return fmt.Sprintf("./testdata/%s-%d.atom", searchTerm, pageNumber)
}

func makeRespFromFile(path string) *http.Response {
	f, err := os.Open(path)
	if err != nil {
		return makeResp(err.Error(), http.StatusInternalServerError)
	}

	okResp := makeResp("200 OK", http.StatusOK)
	okResp.Body = f
	return okResp
}

func resultsPageFromFile(searchTerm string, pageNumber int64) *arxiv.ResultsPage {
	path := resultsPagePath(searchTerm, pageNumber)
	blob, err := ioutil.ReadFile(path)
	if err != nil {
		return nil
	}
	feed := new(Feed)
	if err := xml.Unmarshal(blob, feed); err != nil {
		return nil
	}
	return &arxiv.ResultsPage{
		Feed:       feed,
		PageNumber: pageNumber,
	}
}
