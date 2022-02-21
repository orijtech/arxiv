package arxiv

type Category string

var categoryTable = getLookupTable()

func (c Category) String() string {
	return lookupCategory(c)
}

func lookupCategory(c Category) string {
	v, ok := categoryTable[c]
	if ok {
		return v
	} else {
		return "Unknown Category: " + string(c)
	}
}

const (
	StatisticsApplications                                      Category = "stat.AP"
	StatisticsComputation                                       Category = "stat.CO"
	StatisticsMachineLearning                                   Category = "stat.ML"
	StatisticsMethodology                                       Category = "stat.ME"
	StatisticsTheory                                            Category = "stat.TH"
	QuantitativeBiologyBiomolecules                             Category = "q-bio.BM"
	QuantitativeBiologyCellBehavior                             Category = "q-bio.CB"
	QuantitativeBiologyGenomics                                 Category = "q-bio.GN"
	QuantitativeBiologyMolecularNetworks                        Category = "q-bio.MN"
	QuantitativeBiologyNeuronsAndCognition                      Category = "q-bio.NC"
	QuantitativeBiologyOther                                    Category = "q-bio.OT"
	QuantitativeBiologyPopulationsAndEvolution                  Category = "q-bio.PE"
	QuantitativeBiologyQuantitativeMethods                      Category = "q-bio.QM"
	QuantitativeBiologySubcellularProcesses                     Category = "q-bio.SC"
	QuantitativeBiologyTissuesAndOrgans                         Category = "q-bio.TO"
	CSArchitecture                                              Category = "cs.AR" // https://org/archive/cs
	CSArtificialIntelligence                                    Category = "cs.AI"
	CSComputationAndLanguage                                    Category = "cs.CL"
	ComputationalComplexity                                     Category = "cs.CC"
	ComputationalEngineeringFinanceAndScience                   Category = "cs.CE"
	CSComputationalGeometry                                     Category = "cs.CG"
	CSGameTheory                                                Category = "cs.GT"
	ComputerVisionAndPatternRecognition                         Category = "cs.CV"
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
	SocialAndInformationNetworks                                Category = "cs.SI"
	SystemsAndControl                                           Category = "cs.SY"
	NonLinearSciencesAdaptationAndSelfOrganizingSystemsCategory Category = "nlin.AO"
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

func getLookupTable() map[Category]string {
	categoryMap := map[Category]string{
		StatisticsApplications:                     "StatisticsApplications",
		StatisticsComputation:                      "StatisticsComputation",
		StatisticsMachineLearning:                  "StatisticsMachineLearning",
		StatisticsMethodology:                      "StatisticsMethodology",
		StatisticsTheory:                           "StatisticsTheory",
		QuantitativeBiologyBiomolecules:            "QuantitativeBiologyBiomolecules",
		QuantitativeBiologyCellBehavior:            "QuantitativeBiologyCellBehavior",
		QuantitativeBiologyGenomics:                "QuantitativeBiologyGenomics",
		QuantitativeBiologyMolecularNetworks:       "QuantitativeBiologyMolecularNetworks",
		QuantitativeBiologyNeuronsAndCognition:     "QuantitativeBiologyNeuronsAndCognition",
		QuantitativeBiologyOther:                   "QuantitativeBiologyOther",
		QuantitativeBiologyPopulationsAndEvolution: "QuantitativeBiologyPopulationsAndEvolution",
		QuantitativeBiologyQuantitativeMethods:     "QuantitativeBiologyQuantitativeMethods",
		QuantitativeBiologySubcellularProcesses:    "QuantitativeBiologySubcellularProcesses",
		QuantitativeBiologyTissuesAndOrgans:        "QuantitativeBiologyTissuesAndOrgans",
		CSArchitecture:                             "CSArchitecture",
		CSArtificialIntelligence:                   "CSArtificialIntelligence",
		CSComputationAndLanguage:                   "CSComputationAndLanguage",
		ComputationalComplexity:                    "ComputationalComplexity",
		ComputationalEngineeringFinanceAndScience:  "ComputationalEngineeringFinanceAndScience",
		CSComputationalGeometry:                    "CSComputationalGeometry",
		CSGameTheory:                               "CSGameTheory",
		ComputerVisionAndPatternRecognition:        "ComputerVisionAndPatternRecognition",
		ComputersAndSociety:                        "ComputersAndSociety",
		CryptographyAndSecurity:                    "CryptographyAndSecurity",
		Databases:                                  "Databases",
		DigitalLibraries:                           "DigitalLibraries",
		DiscreteMathematics:                        "DiscreteMathematics",
		DistributedParallelAndClusterComputing:     "DistributedParallelAndClusterComputing",
		CSGeneralLiterature:                        "CSGeneralLiterature",
		CSGraphics:                                 "CSGraphics",
		HumanComputerInteraction:                   "HumanComputerInteraction",
		CSInformationRetrieval:                     "CSInformationRetrieval",
		CSInformationTheory:                        "CSInformationTheory",
		CSLearning:                                 "CSLearning",
		CSLogic:                                    "CSLogic",
		CSMathematicalSoftware:                     "CSMathematicalSoftware",
		MultiagentSystems:                          "MultiagentSystems",
		CSMultimedia:                               "CSMultimedia",
		NetworkAndInternetArchitecture:             "NetworkAndInternetArchitecture",
		NeuralAndEvolutionaryComputing:             "NeuralAndEvolutionaryComputing",
		CSNumericalAnalysis:                        "CSNumericalAnalysis",
		OperatingSystems:                           "OperatingSystems",
		CSOther:                                    "CSOther",
		CSPerformance:                              "CSPerformance",
		ProgrammingLanguages:                       "ProgrammingLanguages",
		CSRobotics:                                 "CSRobotics",
		SoftwareEngineering:                        "SoftwareEngineering",
		CSSound:                                    "CSSound",
		SymbolicComputation:                        "SymbolicComputation",
		SocialAndInformationNetworks:               "SocialAndInformationNetworks",
		SystemsAndControl:                          "SystemsAndControl",
		NonLinearSciencesAdaptationAndSelfOrganizingSystemsCategory: "NonLinearSciencesAdaptationAndSelfOrganizingSystemsCategory",
		NonLinearSciencesCellularAutomataAndLatticeGases:            "NonLinearSciencesCellularAutomataAndLatticeGases",
		NonLinearSciencesChaoticDynamics:                            "NonLinearSciencesChaoticDynamics",
		ExactlySolvableAndIntegrableSytems:                          "ExactlySolvableAndIntegrableSytems",
		PatternFormationAndSolutions:                                "PatternFormationAndSolutions",
		AlgebraicGeometry:                                           "AlgebraicGeometry",
		AlgebraicTopology:                                           "AlgebraicTopology",
		AnalysisOfPDEs:                                              "AnalysisOfPDEs",
		CategoryTheory:                                              "CategoryTheory",
		ClassicalAnalysisAndODEs:                                    "ClassicalAnalysisAndODEs",
		Combinatorics:                                               "Combinatorics",
		CommutativeAlgebra:                                          "CommutativeAlgebra",
		ComplexVariables:                                            "ComplexVariables",
		DifferentialGeometry:                                        "DifferentialGeometry",
		DynamicalSystems:                                            "DynamicalSystems",
		FunctionalAnalysis:                                          "FunctionalAnalysis",
		GeneralMathematics:                                          "GeneralMathematics",
		GeneralTopology:                                             "GeneralTopology",
		GeometricTopology:                                           "GeometricTopology",
		GroupTheory:                                                 "GroupTheory",
		MathsHistoryAndOverview:                                     "MathsHistoryAndOverview",
		MathsInformationTheory:                                      "MathsInformationTheory",
		KTheoryAndHomology:                                          "KTheoryAndHomology",
		MathsLogic:                                                  "MathsLogic",
		MathsMathematicalPhysics:                                    "MathsMathematicalPhysics",
		MetricGeometry:                                              "MetricGeometry",
		NumberTheory:                                                "NumberTheory",
		MathsNumericalAnalysis:                                      "MathsNumericalAnalysis",
		OperatorAlgebras:                                            "OperatorAlgebras",
		MathsOptimizationAndControl:                                 "MathsOptimizationAndControl",
		Probability:                                                 "Probability",
		QuantumAlgebra:                                              "QuantumAlgebra",
		RepresentationTheory:                                        "RepresentationTheory",
		RingsAndAlgebra:                                             "RingsAndAlgebra",
		MathsSpectralTheory:                                         "MathsSpectralTheory",
		MathsStatics:                                                "MathsStatics",
		SymplecticGeometry:                                          "SymplecticGeometry",
		Astrophysics:                                                "Astrophysics",
		PhysicsDisorderedSystemsAndNeuralNetworks:                   "PhysicsDisorderedSystemsAndNeuralNetworks",
		PhysicsMesoscopicSystemsAndQuantumHallEffect:                "PhysicsMesoscopicSystemsAndQuantumHallEffect",
		PhysicsMaterialsScience:                                     "PhysicsMaterialsScience",
		PhysicsOther:                                                "PhysicsOther",
		PhysicsSoftCondensedMatter:                                  "PhysicsSoftCondensedMatter",
		PhysicsStatisticalMechanics:                                 "PhysicsStatisticalMechanics",
		PhysicsStronglyCorrelatedElectrons:                          "PhysicsStronglyCorrelatedElectrons",
		PhysicsSuperconductivity:                                    "PhysicsSuperconductivity",
		GeneralRelativityAndQuantumCosmology:                        "GeneralRelativityAndQuantumCosmology",
		HighEneryPhysicsExperiment:                                  "HighEneryPhysicsExperiment",
		HighEneryPhysicsLattice:                                     "HighEneryPhysicsLattice",
		HighEneryPhysicsPhenomenology:                               "HighEneryPhysicsPhenomenology",
		HighEneryPhysicsTheory:                                      "HighEneryPhysicsTheory",
		MathematicalPhysics:                                         "MathematicalPhysics",
		NuclearExperiment:                                           "NuclearExperiment",
		NuclearTheory:                                               "NuclearTheory",
		AcceleratorPhysics:                                          "AcceleratorPhysics",
		AtmoshpericAndOceanicPhysics:                                "AtmoshpericAndOceanicPhysics",
		AtomicPhysics:                                               "AtomicPhysics",
		AtomicAndMolecularClusters:                                  "AtomicAndMolecularClusters",
		BiologicalPhysics:                                           "BiologicalPhysics",
		ChemicalPhysics:                                             "ChemicalPhysics",
		ClassicalPhysics:                                            "ClassicalPhysics",
		ComputationalPhysics:                                        "ComputationalPhysics",
		DataAnalysisStatisticsAndProbability:                        "DataAnalysisStatisticsAndProbability",
		FluidDynamics:                                               "FluidDynamics",
		GeneralPhysics:                                              "GeneralPhysics",
		Geophysics:                                                  "Geophysics",
		HistoryOfPhysics:                                            "HistoryOfPhysics",
		InstrumentationAndDetectors:                                 "InstrumentationAndDetectors",
		MedicalPhysics:                                              "MedicalPhysics",
		Optics:                                                      "Optics",
		PhysicsEducation:                                            "PhysicsEducation",
		PhysicsAndSociety:                                           "PhysicsAndSociety",
		PlasmaPhysics:                                               "PlasmaPhysics",
		PopularPhysics:                                              "PopularPhysics",
		SpacePhysics:                                                "SpacePhysics",
		QuantumPhysics:                                              "QuantumPhysics",
	}
	return categoryMap
}
