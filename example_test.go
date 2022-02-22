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
	"context"
	"fmt"
	"log"
)

func Example_client_Search_simple() {
	resChan, cancel, err := arxiv.Search(context.Background(), &arxiv.Query{
		Terms:         "deep learning",
		MaxPageNumber: 5,
	})
	if err != nil {
		log.Fatal(err)
	}

	for resPage := range resChan {
		if err := resPage.Err; err != nil {
			fmt.Printf("#%d err: %v", resPage.PageNumber, err)
			continue
		}

		fmt.Printf("#%d\n", resPage.PageNumber)
		feed := resPage.Feed
		fmt.Printf("\tTitle: %s\n\tID: %s\n\tAuthor: %#v\n\tUpdated: %#v\n", feed.Title, feed.ID, feed.Author, feed.Updated)
		for i, entry := range feed.Entry {
			fmt.Printf("\n\t\tEntry: #%d Title: %s ID: %s\n\t\tSummary: %s\n\t\tContent: %#v\n\t\tUpdated: %#v\n\t\tLinks: %#v\n",
				i, entry.Title, entry.ID, entry.Summary.Body, entry.Content, entry.Updated, entry.Link,
			)
		}
		if resPage.PageNumber >= 2 {
			cancel()
		}
	}
}

func Example_client_Search_complex() {
	resChan, cancel, err := arxiv.Search(context.Background(), &arxiv.Query{
		Filters: []*arxiv.Filter{
			{
				Op: arxiv.OpOR,
				Fields: []*arxiv.Field{
					{Title: "architecture"},
					{Category: arxiv.CSGameTheory},
					{Comment: "biological"},
				},
			},
		},
		MaxPageNumber: 2,
	})
	if err != nil {
		log.Fatal(err)
	}

	for resPage := range resChan {
		if err := resPage.Err; err != nil {
			fmt.Printf("#%d err: %v", resPage.PageNumber, err)
			continue
		}

		fmt.Printf("#%d\n", resPage.PageNumber)
		feed := resPage.Feed
		fmt.Printf("\tTitle: %s\n\tID: %s\n\tAuthor: %#v\n\tUpdated: %#v\n", feed.Title, feed.ID, feed.Author, feed.Updated)
		for i, entry := range feed.Entry {
			fmt.Printf("\n\t\tEntry: #%d Title: %s ID: %s\n\t\tSummary: %s\n\t\tContent: %#v\n\t\tUpdated: %#v\n\t\tLinks: %#v\n",
				i, entry.Title, entry.ID, entry.Summary.Body, entry.Content, entry.Updated, entry.Link,
			)
		}
		if resPage.PageNumber >= 2 {
			cancel()
		}
	}
}
