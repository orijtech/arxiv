# arxiv

Patched version of the stalled go-lang client. 
* Merged bugfix for id_list parameter
* Refactored some code. 


Go API client for arxiv.org. It supports simple as well as advanced searches with filters.


## Usage
Sample usage can be found in file [example_test.go](./example_test.go)
Or see below:
* Preamble:
```go
package main

import (
    "fmt"
    "log"

    "github.com/marvin-hansen/arxiv/v1"
)
```

## Examples
* Simple search
```go
func simpleSearch() {
	resChan, cancel, err := arxiv.Search(&arxiv.Query{
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
```

* Advanced/complex search
```go
func advancedSearch() {
	resChan, cancel, err := arxiv.Search(&arxiv.Query{
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
```
