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

// Category

type Category string

const (
	StatisticsApplications                     Category = "stat.AP"
	StatisticsComputation                      Category = "stat.CO"
	StatisticsMachineLearning                  Category = "stat.ML"
	StatisticsMethodology                      Category = "stat.ME"
	StatisticsTheory                           Category = "stat.TH"
	QuantitativeBiologyBiomolecules            Category = "q-bio.BM"
	QuantitativeBiologyCellBehavior            Category = "q-bio.CB"
	QuantitativeBiologyGenomics                Category = "q-bio.GN"
	QuantitativeBiologyMolecularNetworks       Category = "q-bio.MN"
	QuantitativeBiologyNeuronsAndCognition     Category = "q-bio.NC"
	QuantitativeBiologyOther                   Category = "q-bio.OT"
	QuantitativeBiologyPopulationsAndEvolution Category = "q-bio.PE"
	QuantitativeBiologyQuantitativeMethods     Category = "q-bio.QM"
	QuantitativeBiologySubcellularProcesses    Category = "q-bio.SC"
	QuantitativeBiologyTissuesAndOrgans        Category = "q-bio.TO"
	CSArchitecture                             Category = "cs.AR"
	CSArtificialIntelligence                   Category = "cs.AI"
	CSComputationAndLanguage                   Category = "cs.CL"
	ComputationalComplexity                    Category = "cs.CC"
	ComputationalEngineeringFinanceAndScience  Category = "cs.CE"
	CSComputationalGeometry                    Category = "cs.CG"
	CSGameTheory                               Category = "cs.GT"
	ComputerVisionAndPatternRecognition        Category = "cs.CV"

	ComputersAndSociety                                         Category = "cs.CY"
	CryptographyAndSecurity                                     Category = "cs.CR"
	Databases                                                   Category = "cs.DB"
	DigitalLibraries                                            Category = "cs.DL"
	DiscreteMathematics                                         Category = "cs.DM"
	DistributedParallelAndClusterComputing                      Category = "cs.DC"
	CSGeneralLiterature                                         Category = "cs.GL"
	CSGraphics                                                  Category = "cs.GR"
	HumanComputerInteraction                                    Category = "cs.HC"
	CSInformationRetrieval                                      Category = "cs.IR"
	CSInformationTheory                                         Category = "cs.IT"
	CSLearning                                                  Category = "cs.LG"
	CSLogic                                                     Category = "cs.LO"
	CSMathematicalSoftware                                      Category = "cs.MS"
	MultiagentSystems                                           Category = "cs.MA"
	CSMultimedia                                                Category = "cs.MM"
	NetworkAndInternetArchitecture                              Category = "cs.NI"
	NeuralAndEvolutionaryComputing                              Category = "cs.NE"
	CSNumericalAnalysis                                         Category = "cs.NA"
	OperatingSystems                                            Category = "cs.OS"
	CSOther                                                     Category = "cs.OH"
	CSPerformance                                               Category = "cs.PF"
	ProgrammingLanguages                                        Category = "cs.PL"
	CSRobotics                                                  Category = "cs.RO"
	SoftwareEngineering                                         Category = "cs.SE"
	CSSound                                                     Category = "cs.SD"
	SymbolicComputation                                         Category = "cs.SC"
	NonLinearSciencesAdaptationAndSelfOrganizingSystemsCategory          = "nlin.AO"
	NonLinearSciencesCellularAutomataAndLatticeGases            Category = "nlin.CG"
	NonLinearSciencesChaoticDynamics                            Category = "nlin.CD"
	ExactlySolvableAndIntegrableSytems                          Category = "nlin.SI"
	PatternFormationAndSolutions                                Category = "nlin.PS"
	AlgebraicGeometry                                           Category = "math.AG"
	AlgebraicTopology                                           Category = "math.AT"
	AnalysisOfPDEs                                              Category = "math.AP"
	CategoryTheory                                              Category = "math.CT"
	ClassicalAnalysisAndODEs                                    Category = "math.CA"
	Combinatorics                                               Category = "math.CO"
	CommutativeAlgebra                                          Category = "math.AC"
	ComplexVariables                                            Category = "math.CV"
	DifferentialGeometry                                        Category = "math.DG"
	DynamicalSystems                                            Category = "math.DS"
	FunctionalAnalysis                                          Category = "math.FA"
	GeneralMathematics                                          Category = "math.GM"
	GeneralTopology                                             Category = "math.GN"
	GeometricTopology                                           Category = "math.GT"
	GroupTheory                                                 Category = "math.GR"
	MathsHistoryAndOverview                                     Category = "math.HO"
	MathsInformationTheory                                      Category = "math.IT"
	KTheoryAndHomology                                          Category = "math.KT"
	MathsLogic                                                  Category = "math.LO"
	MathsMathematicalPhysics                                    Category = "math.MP"
	MetricGeometry                                              Category = "math.MG"
	NumberTheory                                                Category = "math.NT"
	MathsNumericalAnalysis                                      Category = "math.NA"
	OperatorAlgebras                                            Category = "math.OA"
	MathsOptimizationAndControl                                 Category = "math.OC"
	Probability                                                 Category = "math.PR"
	QuantumAlgebra                                              Category = "math.QA"
	RepresentationTheory                                        Category = "math.RT"
	RingsAndAlgebra                                             Category = "math.RA"
	MathsSpectralTheory                                         Category = "math.SP"
	MathsStatics                                                Category = "math.ST"
	SymplecticGeometry                                          Category = "math.SG"
	Astrophysics                                                Category = "astro-ph"
	PhysicsDisorderedSystemsAndNeuralNetworks                   Category = "cond-mat.dis-nn"
	PhysicsMesoscopicSystemsAndQuantumHallEffect                Category = "cond-mat.mes-hall"
	PhysicsMaterialsScience                                     Category = "cond-mat.mtrl-sci"
	PhysicsOther                                                Category = "cond-mat.other"
	PhysicsSoftCondensedMatter                                  Category = "cond-mat.soft"
	PhysicsStatisticalMechanics                                 Category = "cond-mat.stat-mech"
	PhysicsStronglyCorrelatedElectrons                          Category = "cond-mat.str-el"
	PhysicsSuperconductivity                                    Category = "cond-mat.sup-con"
	GeneralRelativityAndQuantumCosmology                        Category = "gr-qc"
	HighEneryPhysicsExperiment                                  Category = "hep-ex"
	HighEneryPhysicsLattice                                     Category = "hep-lat"
	HighEneryPhysicsPhenomenology                               Category = "hep-ph"
	HighEneryPhysicsTheory                                      Category = "hep-th"
	MathematicalPhysics                                         Category = "math-ph"
	NuclearExperiment                                           Category = "nucl-ex"
	NuclearTheory                                               Category = "nucl-th"
	AcceleratorPhysics                                          Category = "physics.acc-ph"
	AtmoshpericAndOceanicPhysics                                Category = "physics.ao-ph"
	AtomicPhysics                                               Category = "physics.atom-ph"
	AtomicAndMolecularClusters                                  Category = "physics.atm-clus"
	BiologicalPhysics                                           Category = "physics.bio-ph"
	ChemicalPhysics                                             Category = "physics.chem-ph"
	ClassicalPhysics                                            Category = "physics.class-ph"
	ComputationalPhysics                                        Category = "physics.comp-ph"
	DataAnalysisStatisticsAndProbability                        Category = "physics.data-an"
	FluidDynamics                                               Category = "physics.flu-dyn"
	GeneralPhysics                                              Category = "physics.gen-ph"
	Geophysics                                                  Category = "physics.geo-ph"
	HistoryOfPhysics                                            Category = "physics.hist-ph"
	InstrumentationAndDetectors                                 Category = "physics.ins-det"
	MedicalPhysics                                              Category = "physics.med-ph"
	Optics                                                      Category = "physics.optics"
	PhysicsEducation                                            Category = "physics.ed-ph"
	PhysicsAndSociety                                           Category = "physics.soc-ph"
	PlasmaPhysics                                               Category = "physiscs.plasm-ph"
	PopularPhysics                                              Category = "physics.pop-ph"
	SpacePhysics                                                Category = "physics.space-ph"
	QuantumPhysics                                              Category = "quant-ph"
)
