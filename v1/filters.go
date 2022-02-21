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
	"fmt"
	"reflect"
	"strings"
)

type Field struct {
	Title            string   `json:"ti,omitempty"`
	Author           string   `json:"au,omitempty"`
	Category         Category `json:"cat,omitempty"`
	Abstract         string   `json:"abs,omitempty"`
	Comment          string   `json:"co,omitempty"`
	ReportNumber     int64    `json:"rn,omitempty"`
	JournalReference string   `json:"jn,omitempty"`

	All string `json:"all,omitempty"`
}

type Operator string

const (
	OpAnd    Operator = "AND"
	OpOR     Operator = "OR"
	OpAndNot Operator = "ANDNOT"
	OpAll    Operator = "all"
)

type Filter struct {
	Op     Operator `json:"op,omitempty"`
	Fields []*Field `json:"fields,omitempty"`
}

var _ = []Filter{
	{
		Op: OpOR,
		Fields: []*Field{
			{Title: "Richard Feynmann"},
			{Author: "Dmitry Vyukov", Category: ""},
			{Comment: "Michael Nielsen", Category: StatisticsMachineLearning},
			{Category: NuclearExperiment},
		},
	},
}

type Filters []*Filter

func (fl Filters) String() string {
	var strs []string
	for _, filter := range fl {
		str := filter.Pairfy()
		if str != "" {
			strs = append(strs, fmt.Sprintf("(%s)", str))
		}
	}
	return strings.Join(strs, fmt.Sprintf(" %s ", OpAnd))
}

// Transform:
// var _ = []Filter{
// 	{
// 		Op: OpAnd,
// 		Fields: []*Field{
// 			{Author: "Adrian DelMaestro"},
// 			{Title: "Checkerboard"},
// 		},
// 	},
// }
// Into:
// au:Adrian%20DelMaestro AND ti:checkerboard
func (f *Filter) Pairfy() string {
	fields := f.Fields
	var pairs []string
	for _, field := range fields {
		val := reflect.ValueOf(field)
		if val.Kind() == reflect.Ptr {
			val = reflect.Indirect(val)
		}
		typ := val.Type()

		var innerFieldOps []string
		nfields := val.NumField()
		for j := 0; j < nfields; j++ {
			jthField := val.Field(j)
			if jthField.Kind() == reflect.Ptr {
				jthField = reflect.Indirect(jthField)
			}

			if jthField.Kind() == reflect.Invalid {
				continue
			}

			jthTypeField := typ.Field(j)
			tag, omitempty, ignore := jsonTag(jthTypeField)
			if ignore {
				continue
			}

			zeroValue := reflect.Zero(jthTypeField.Type)
			storedValue := jthField.Interface()

			// FYI custom blank checks using untyped zero values e.g
			// "", 0, false
			// fail for type aliases and int64, int32, uint, uint32 etc.
			// See https://play.golang.org/p/ID8L4kVO10.
			isZero := storedValue == zeroValue.Interface()
			if omitempty && isZero {
				continue
			}

			innerFieldOps = append(innerFieldOps, fmt.Sprintf("%s:%v", tag, storedValue))
		}

		if len(innerFieldOps) > 0 {
			// Within the same field, join operators using the AND operator
			// e.g {Title: "Foo", Comment; "Golang"} ==> "ti:Foo+AND+co:Golang"
			finalOp := strings.Join(innerFieldOps, " AND ")
			if len(innerFieldOps) > 1 {
				finalOp = "(" + finalOp + ")"
			}
			pairs = append(pairs, finalOp)
		}
	}

	return strings.Join(pairs, fmt.Sprintf(" %s ", f.Op))
}

func jsonTag(v reflect.StructField) (tag string, omitempty, ignore bool) {
	jsonTag := v.Tag.Get("json")
	if jsonTag == "" {
		return v.Name, false, false
	}

	splits := strings.Split(jsonTag, ",")
	if len(splits) == 0 {
		return "", false, false
	}
	tag, instrs := splits[0], splits[1:]
	instrIndex := make(map[string]bool)
	for _, instr := range instrs {
		instrIndex[instr] = true
	}

	_, omitempty = instrIndex["omitempty"]
	_, ignore = instrIndex["-"]
	return tag, omitempty, ignore
}
