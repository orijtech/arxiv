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
	"testing"
)

func TestFilters(t *testing.T) {
	tests := [...]struct {
		filter *arxiv.Filter
		want   string
	}{
		0: {
			want: "au:Adrian DelMaestro AND ti:Checkerboard",
			filter: &arxiv.Filter{
				Op: arxiv.OpAnd,
				Fields: []*arxiv.Field{
					{Author: "Adrian DelMaestro"},
					{Title: "Checkerboard"},
				},
			},
		},

		1: {
			want: "(au:Richard Mueller AND co:Berkeley) OR ti:Nuclear OR cat:cs.GT",
			filter: &arxiv.Filter{
				Op: arxiv.OpOR,
				Fields: []*arxiv.Field{
					{Author: "Richard Mueller", Comment: "Berkeley"},
					{Title: "Nuclear"},
					{Category: arxiv.CSGameTheory},
				},
			},
		},
	}

	for i, tt := range tests {
		got, want := tt.filter.Pairfy(), tt.want
		if got != want {
			t.Errorf("#%d\ngot: %q\nwant:%q\n", i, got, want)
		}
	}
}
