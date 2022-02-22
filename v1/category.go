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

// Category Taxonomy. Updated Feb/22 2022
// https://arxiv.org/category_taxonomy
const ( // Computer Science
	ArtificialIntelligence                    Category = "cs.AI" // Covers all areas of AI except Vision, Robotics, Machine Learning, Multiagent Systems, and Computation and Language (Natural Language Processing), which have separate subject areas.
	HardwareArchitecture                      Category = "cs.AR" // Covers systems organization and hardware architecture. Roughly includes material in ACM Subject Classes C.0, C.1, and C.5.
	ComputationalComplexity                   Category = "cs.CC" // Covers models of computation, complexity classes, structural complexity, complexity tradeoffs, upper and lower bounds.
	ComputationalEngineeringFinanceAndScience Category = "cs.CE" // Covers applications of computer science to the mathematical modeling of complex systems in the fields of science, engineering, and finance.
	ComputationalGeometry                     Category = "cs.CG" // Roughly includes material in ACM Subject Classes I.3.5 and F.2.2.
	ComputationAndLanguage                    Category = "cs.CL" // Covers natural language processing. Roughly includes material in ACM Subject Class I.2.7.
	CryptographyAndSecurity                   Category = "cs.CR" // Covers all areas of cryptography and security including authentication, public key cryptosytems, proof-carrying code, etc. Roughly includes material in ACM Subject Classes D.4.6 and E.3.
	ComputerVisionAndPatternRecognition       Category = "cs.CV" // Covers image processing, computer vision, pattern recognition, and scene understanding. Roughly includes material in ACM Subject Classes I.2.10, I.4, and I.5.
	ComputersAndSociety                       Category = "cs.CY" // Covers impact of computers on society, computer ethics, information technology and public policy, legal aspects of computing, computers and education. Roughly includes material in ACM Subject Classes K.0, K.2, K.3, K.4, K.5, and K.7.
	Databases                                 Category = "cs.DB" // Covers database management, datamining, and data processing. Roughly includes material in ACM Subject Classes E.2, E.5, H.0, H.2, and J.1.
	DistributedParallelAndClusterComputing    Category = "cs.DC" // Covers fault-tolerance, distributed algorithms, stabilility, parallel computation, and cluster computing. Roughly includes material in ACM Subject Classes C.1.2, C.1.4, C.2.4, D.1.3, D.4.5, D.4.7, E.1.
	DigitalLibraries                          Category = "cs.DL" // Covers all aspects of the digital library design and document and text creation. Note that there will be some overlap with Information Retrieval (which is a separate subject area). Roughly includes material in ACM Subject Classes H.3.5, H.3.6, H.3.7, I.7.
	DiscreteMathematics                       Category = "cs.DM" // Covers combinatorics, graph theory, applications of probability. Roughly includes material in ACM Subject Classes G.2 and G.3.
	DataStructuresAndAlgorithms               Category = "cs.DS" // Covers data structures and analysis of algorithms. Roughly includes material in ACM Subject Classes E.1, E.2, F.2.1, and F.2.2.
	EmergingTechnologies                      Category = "cs.ET" // Covers approaches to information processing (computing, communication, sensing) and bio-chemical analysis based on alternatives to silicon CMOS-based technologies, such as nanoscale electronic, photonic, spin-based, superconducting, mechanical, bio-chemical and quantum technologies (this list is not exclusive)
	FormalLanguagesAndAutomataTheory          Category = "cs.FL" // Covers automata theory, formal language theory, grammars, and combinatorics on words. This roughly corresponds to ACM Subject Classes F.1.1, and F.4.3. Papers dealing with computational complexity should go to cs.CC; papers dealing with logic should go to cs.LO.
	GeneralLiterature                         Category = "cs.GL" // Covers introductory material, survey material, predictions of future trends, biographies, and miscellaneous computer-science related material. Roughly includes all of ACM Subject Class A, except it does not include conference proceedings (which will be listed in the appropriate subject area).
	Graphics                                  Category = "cs.GR" // Covers all aspects of computer graphics. Roughly includes material in all of ACM Subject Class I.3, except that I.3.5 is is likely to have Computational Geometry as the primary subject area.
	ComputerScienceAndGameTheory              Category = "cs.GT" // Covers all theoretical and applied aspects at the intersection of computer science and game theory, including work in mechanism design, learning in games (which may overlap with Learning), foundations of agent modeling in games (which may overlap with Multiagent systems), coordination, specification and formal methods for non-cooperative computational environments. T
	HumanComputerInteraction                  Category = "cs.HC" // Covers human factors, user interfaces, and collaborative computing. Roughly includes material in ACM Subject Classes H.1.2 and all of H.5, except for H.5.1, which is more likely to have Multimedia as the primary subject area.
	InformationRetrieval                      Category = "cs.IR" // Covers indexing, dictionaries, retrieval, content and analysis. Roughly includes material in ACM Subject Classes H.3.0, H.3.1, H.3.2, H.3.3, and H.3.4.
	InformationTheory                         Category = "cs.IT" // Covers theoretical and experimental aspects of information theory and coding. Includes material in ACM Subject Class E.4 and intersects with H.1.1.
	MachineLearning                           Category = "cs.LG" // Papers on all aspects of machine learning research (supervised, unsupervised, reinforcement learning, bandit problems, and so on) including also robustness, explanation, fairness, and methodology. cs.LG is also an appropriate primary category for applications of machine learning methods.
	LogicInComputerScience                    Category = "cs.LO" // Covers all aspects of logic in computer science, including finite model theory, logics of programs, modal logic, and program verification. Programming language semantics should have Programming Languages as the primary subject area.
	MultiagentSystems                         Category = "cs.MA" // Covers multiagent systems, distributed artificial intelligence, intelligent agents, coordinated interactions. and practical applications. Roughly covers ACM Subject Class I.2.11.
	Multimedia                                Category = "cs.MM" // Roughly includes material in ACM Subject Class H.5.1.
	MathematicalSoftware                      Category = "cs.MS" // Roughly includes material in ACM Subject Class G.4.
	NumericalAnalysis                         Category = "cs.NA" // cs.NA is an alias for math.NA. Roughly includes material in ACM Subject Class G.1.
	NeuralAndEvolutionaryComputing            Category = "cs.NE" // Covers neural networks, connectionism, genetic algorithms, artificial life, adaptive behavior. Roughly includes some material in ACM Subject Class C.1.3, I.2.6, I.5.
	NetworkAndInternetArchitecture            Category = "cs.NI" // Covers all aspects of computer communication networks, including network architecture and design, network protocols, and internetwork standards (like TCP/IP). Also includes topics, such as web caching, that are directly relevant to Internet architecture and performance.
	CSOther                                   Category = "cs.OH" // This is the classification to use for documents that do not fit anywhere else.
	OperatingSystems                          Category = "cs.OS" // Roughly includes material in ACM Subject Classes D.4.1, D.4.2., D.4.3, D.4.4, D.4.5, D.4.7, and D.4.9.
	Performance                               Category = "cs.PF" // Covers performance measurement and evaluation, queueing, and simulation. Roughly includes material in ACM Subject Classes D.4.8 and K.6.2.
	ProgrammingLanguages                      Category = "cs.PL" // Covers programming language semantics, language features, programming approaches (such as object-oriented programming, functional programming, logic programming). Also includes material on compilers oriented towards programming languages; other material on compilers may be more appropriate in Architecture (AR). Roughly includes material in ACM Subject Classes D.1 and D.3.
	Robotics                                  Category = "cs.RO" // Roughly includes material in ACM Subject Class I.2.9.
	SymbolicComputation                       Category = "cs.SC" // Roughly includes material in ACM Subject Class I.1.
	Sound                                     Category = "cs.SD" // Covers all aspects of computing with sound, and sound as an information channel. Includes models of sound, analysis and synthesis, audio user interfaces, sonification of data, computer music, and sound signal processing.
	SoftwareEngineering                       Category = "cs.SE" // Covers design tools, software metrics, testing and debugging, programming environments, etc. Roughly includes material in all of ACM Subject Classes D.2, except that D.2.4 (program verification) should probably have Logics in Computer Science as the primary subject area.
	SocialAndInformationNetworks              Category = "cs.SI" // Covers the design, analysis, and modeling of social and information networks, including their applications for on-line information access, communication, and interaction, and their roles as datasets in the exploration of questions in these and other domains, including connections to the social and biological sciences.
	CSSystemsAndControl                       Category = "cs.SY" //  cs.SY is an alias for eess.SY. This section includes theoretical and experimental research covering all facets of automatic control systems. The section is focused on methods of control system analysis and design using tools of modeling, simulation and optimization. Specific areas of research include nonlinear, distributed, adaptive, stochastic and robust control in addition to hybrid and discrete event systems.
) // Computer Science

const ( // Economics
	Econometrics         Category = "econ.EM" // Econometric Theory, Micro-Econometrics, Macro-Econometrics, Empirical Content of Economic Relations discovered via New Methods, Methodological Aspects of the Application of Statistical Inference to Economic Data.
	GeneralEconomics     Category = "econ.GN" // General methodological, applied, and empirical contributions to economics.
	TheoreticalEconomics Category = "econ.TH" // Includes theoretical contributions to Contract Theory, Decision Theory, Game Theory, General Equilibrium, Growth, Learning and Evolution, Macroeconomics, Market and Mechanism Design, and Social Choice.
) // Economics

const ( // Electrical Engineering and Systems Science
	AudioSndSpeechProcessing Category = "eess.AS" // Theory and methods for processing signals representing audio, speech, and language, and their applications. T
	ImageAndVideoProcessing  Category = "eess.IV" // Theory, algorithms, and architectures for the formation, capture, processing, communication, analysis, and display of images, video, and multidimensional signals in a wide variety of applications.
	SignalProcessing         Category = "eess.SP" // Theory, algorithms, performance analysis and applications of signal and data analysis, including physical modeling, processing, detection and parameter estimation, learning, mining, retrieval, and information extraction. The term "signal" includes speech, audio, sonar, radar, geophysical, physiological, (bio-) medical, image, video, and multimodal natural and man-made signals, including communication signals and data.
	SystemsAndControl        Category = "eess.SY" // This section includes theoretical and experimental research covering all facets of automatic control systems. The section is focused on methods of control system analysis and design using tools of modeling, simulation and optimization.
) // Electrical Engineering and Systems Science

const ( // Mathematics
	CommutativeAlgebra          Category = "math.AC" // Commutative rings, modules, ideals, homological algebra, computational aspects, invariant theory, connections to algebraic geometry and combinatorics
	AlgebraicGeometry           Category = "math.AG" // Algebraic varieties, stacks, sheaves, schemes, moduli spaces, complex geometry, quantum cohomology
	AnalysisOfPDEs              Category = "math.AP" // Existence and uniqueness, boundary conditions, linear and non-linear operators, stability, soliton theory, integrable PDE's, conservation laws, qualitative dynamics
	AlgebraicTopology           Category = "math.AT" // Homotopy theory, homological algebra, algebraic treatments of manifolds
	ClassicalAnalysisAndODEs    Category = "math.CA" // Special functions, orthogonal polynomials, harmonic analysis, ODE's, differential relations, calculus of variations, approximations, expansions, asymptotics
	Combinatorics               Category = "math.CO" // Discrete mathematics, graph theory, enumeration, combinatorial optimization, Ramsey theory, combinatorial game theory
	CategoryTheory              Category = "math.CT" // Enriched categories, topoi, abelian categories, monoidal categories, homological algebra
	ComplexVariables            Category = "math.CV" // Holomorphic functions, automorphic group actions and forms, pseudoconvexity, complex geometry, analytic spaces, analytic sheaves
	DifferentialGeometry        Category = "math.DG" // Complex, contact, Riemannian, pseudo-Riemannian and Finsler geometry, relativity, gauge theory, global analysis
	DynamicalSystems            Category = "math.DS" // Dynamics of differential equations and flows, mechanics, classical few-body problems, iterations, complex dynamics, delayed differential equations
	FunctionalAnalysis          Category = "math.FA" // Banach spaces, function spaces, real functions, integral transforms, theory of distributions, measure theory
	GeneralMathematics          Category = "math.GM" // Mathematical material of general interest, topics not covered elsewhere
	GeneralTopology             Category = "math.GN" // Continuum theory, point-set topology, spaces with algebraic structure, foundations, dimension theory, local and global properties
	GeometricTopology           Category = "math.GT" // Finite groups, topological groups, representation theory, cohomology, classification and structure
	GroupTheory                 Category = "math.GR" // Manifolds, orbifolds, polyhedra, cell complexes, foliations, geometric structures
	MathsHistoryAndOverview     Category = "math.HO" // Biographies, philosophy of mathematics, mathematics education, recreational mathematics, communication of mathematics, ethics in mathematics
	MathsInformationTheory      Category = "math.IT" // math.IT is an alias for cs.IT. Covers theoretical and experimental aspects of information theory and coding.
	KTheoryAndHomology          Category = "math.KT" // Algebraic and topological K-theory, relations with topology, commutative algebra, and operator algebras
	MathsLogic                  Category = "math.LO" // Logic, set theory, point-set topology, formal mathematics
	MathsMathematicalPhysics    Category = "math.MP" // Euclidean, hyperbolic, discrete, convex, coarse geometry, comparisons in Riemannian geometry, symmetric spaces
	MetricGeometry              Category = "math.MG" // math.MP is an alias for math-ph. Articles in this category focus on areas of research that illustrate the application of mathematics to problems in physics, develop mathematical methods for such applications, or provide mathematically rigorous formulations of existing physical theories.
	MathsNumericalAnalysis      Category = "math.NA" // Numerical algorithms for problems in analysis and algebra, scientific computation
	NumberTheory                Category = "math.NT" // Prime numbers, diophantine equations, analytic number theory, algebraic number theory, arithmetic geometry, Galois theory
	OperatorAlgebras            Category = "math.OA" // Algebras of operators on Hilbert space, C^*-algebras, von Neumann algebras, non-commutative geometry
	MathsOptimizationAndControl Category = "math.OC" // Operations research, linear programming, control theory, systems theory, optimal control, game theory
	Probability                 Category = "math.PR" // Theory and applications of probability and stochastic processes: e.g. central limit theorems, large deviations, stochastic differential equations, models from statistical mechanics, queuing theory
	QuantumAlgebra              Category = "math.QA" // Quantum groups, skein theories, operadic and diagrammatic algebra, quantum field theory
	RingsAndAlgebra             Category = "math.RA" // Non-commutative rings and algebras, non-associative algebras, universal algebra and lattice theory, linear algebra, semigroups
	RepresentationTheory        Category = "math.RT" // Linear representations of algebras and groups, Lie theory, associative algebras, multilinear algebra
	SymplecticGeometry          Category = "math.SG" // Hamiltonian systems, symplectic flows, classical integrable systems
	MathsSpectralTheory         Category = "math.SP" // Schrodinger operators, operators on manifolds, general differential operators, numerical studies, integral operators, discrete models, resonances, non-self-adjoint operators, random operators/matrices
	MathsStatics                Category = "math.ST" // Applied, computational and theoretical statistics: e.g. statistical inference, regression, time series, multivariate analysis, data analysis, Markov chain Monte Carlo, design of experiments, case studies
) // Mathematics

const ( // Physics
	Physics                                                     Category = "physics"
	Astrophysics                                                Category = "astro-ph"
	CosmologyAndNongalacticAstrophysics                         Category = "astro-ph.CO"        // Phenomenology of early universe, cosmic microwave background, cosmological parameters, primordial element abundances, extragalactic distance scale, large-scale structure of the universe. Groups, superclusters, voids, intergalactic medium. P
	EarthAndPlanetaryAstrophysics                               Category = "astro-ph.EP"        // Interplanetary medium, planetary physics, planetary astrobiology, extrasolar planets, comets, asteroids, meteorites. Structure and formation of the solar system
	AstrophysicsOfGalaxies                                      Category = "astro-ph.GA"        // Phenomena pertaining to galaxies or the Milky Way. Star clusters, HII regions and planetary nebulae, the interstellar medium, atomic and molecular clouds, dust. Stellar populations.
	HighEnergyAstrophysicalPhenomena                            Category = "astro-ph.HE"        // Cosmic ray production, acceleration, propagation, detection. Gamma ray astronomy and bursts, X-rays, charged particles, supernovae and other explosive phenomena, stellar remnants and accretion systems, jets, microquasars, neutron stars, pulsars, black holes
	InstrumentationAndMethodsForAstrophysics                    Category = "astro-ph.IM"        // Detector and telescope design, experiment proposals. Laboratory Astrophysics. Methods for data analysis, statistical methods. Software, database design
	SolarAndStellarAstrophysics                                 Category = "astro-ph.SR"        // White dwarfs, brown dwarfs, cataclysmic variables. Star formation and protostellar systems, stellar astrobiology, binary and multiple systems of stars, stellar evolution and structure, coronas. Central stars of planetary nebulae. Helioseismology, solar neutrinos, production and detection of gravitational radiation from stellar systems
	CondensedMatter                                             Category = "cond-mat"           // CondensedMatter
	DisorderedSystemsAndNeuralNetworks                          Category = "cond-mat.dis-nn"    // Glasses and spin glasses; properties of random, aperiodic and quasiperiodic systems; transport in disordered media; localization; phenomena mediated by defects and disorder; neural networks
	MesoscaleAndNanoscalePhysics                                Category = "cond-mat.mes-hall"  // Semiconducting nanostructures: quantum dots, wires, and wells. Single electronics, spintronics, 2d electron gases, quantum Hall effect, nanotubes, graphene, plasmonic nanostructures
	MaterialsScience                                            Category = "cond-mat.mtrl-sci"  // Techniques, synthesis, characterization, structure. Structural phase transitions, mechanical properties, phonons. Defects, adsorbates, interfaces
	PhysicsOther                                                Category = "cond-mat.other"     // Work in condensed matter that does not fit into the other cond-mat classifications
	QuantumGases                                                Category = "cond-mat.quant-gas" // Ultracold atomic and molecular gases, Bose-Einstein condensation, Feshbach resonances, spinor condensates, optical lattices, quantum simulation with cold atoms and molecules, macroscopic interference phenomena
	SoftCondensedMatter                                         Category = "cond-mat.soft"      //Membranes, polymers, liquid crystals, glasses, colloids, granular matter
	StatisticalMechanics                                        Category = "cond-mat.stat-mech" // Phase transitions, thermodynamics, field theory, non-equilibrium phenomena, renormalization group and scaling, integrable models, turbulence
	StronglyCorrelatedElectrons                                 Category = "cond-mat.str-el"    // Quantum magnetism, non-Fermi liquids, spin liquids, quantum criticality, charge density waves, metal-insulator transitions
	Superconductivity                                           Category = "cond-mat.sup-con"   // Superconductivity: theory, models, experiment. Superflow in helium
	GeneralRelativityAndQuantumCosmology                        Category = "gr-qc"              // General Relativity and Quantum Cosmology Areas of gravitational physics, including experiments and observations related to the detection and interpretation of gravitational waves, experimental tests of gravitational theories, computational general relativity, relativistic astrophysics, solutions to Einstein's equations and their properties, alternative theories of gravity, classical and quantum cosmology, and quantum gravity.
	HighEnergyPhysic                                            Category = "hep"                // High Energy Physics
	HighEnergyPhysicsExperiment                                 Category = "hep-ex"             // High Energy Physics - Experiment
	HighEnergyPhysicsLattice                                    Category = "hep-lat"            // Lattice field theory. Phenomenology from lattice field theory. Algorithms for lattice field theory. Hardware for lattice field theory.
	HighEnergyPhysicsPhenomenology                              Category = "hep-ph"             // Theoretical particle physics and its interrelation with experiment. Prediction of particle physics observables: models, effective field theories, calculation techniques. Particle physics: analysis of theory through experimental results.
	HighEnergyPhysicsTheory                                     Category = "hep-th"             // Formal aspects of quantum field theory. String theory, supersymmetry and supergravity.
	MathematicalPhysics                                         Category = "math-ph"            // Articles in this category focus on areas of research that illustrate the application of mathematics to problems in physics, develop mathematical methods for such applications, or provide mathematically rigorous formulations of existing physical theories.
	NonLinearSciences                                           Category = "nlin"               // NonLinearSciences
	NonLinearSciencesAdaptationAndSelfOrganizingSystemsCategory Category = "nlin.AO"            // Adaptation, self-organizing systems, statistical physics, fluctuating systems, stochastic processes, interacting particle systems, machine learning
	NonLinearSciencesChaoticDynamics                            Category = "nlin.CD"            // Dynamical systems, chaos, quantum chaos, topological dynamics, cycle expansions, turbulence, propagation
	NonLinearSciencesCellularAutomataAndLattice                 Category = "nlin.CG"            //Computational methods, time series analysis, signal processing, wavelets, lattice gases
	PatternFormationAndSolutions                                Category = "nlin.PS"            // Pattern formation, coherent structures, solitons
	ExactlySolvableAndIntegrableSystems                         Category = "nlin.SI"            // Exactly solvable systems, integrable PDEs, integrable ODEs, Painleve analysis, integrable discrete maps, solvable lattice models, integrable quantum systems
	NuclearExperiment                                           Category = "nucl-ex"            // Nuclear Experiment Results from experimental nuclear physics including the areas of fundamental interactions, measurements at low- and medium-energy, as well as relativistic heavy-ion collisions.
	NuclearTheory                                               Category = "nucl-th"            // Nuclear Theory Theory of nuclear structure covering wide area from models of hadron structure to neutron stars. Nuclear equation of states at different external conditions.
	AcceleratorPhysics                                          Category = "physics.acc-ph"     // Accelerator theory and simulation. Accelerator technology. Accelerator experiments. Beam Physics. Accelerator design and optimization. Advanced accelerator concepts. Radiation sources including synchrotron light sources and free electron lasers. Applications of accelerators.
	AtmoshpericAndOceanicPhysics                                Category = "physics.ao-ph"      // Atmospheric and oceanic physics and physical chemistry, biogeophysics, and climate science
	AppliedPhysics                                              Category = "physics.app-ph"     // Applications of physics to new technology, including electronic devices, optics, photonics, microwaves, spintronics, advanced materials, metamaterials, nanotechnology, and energy sciences.
	AtomicAndMolecularClusters                                  Category = "physics.atm-clus"   // Atomic and molecular clusters, nanoparticles: geometric, electronic, optical, chemical, magnetic properties, shell structure, phase transitions, optical spectroscopy, mass spectrometry, photoelectron spectroscopy, ionization potential, electron affinity, interaction with intense light pulses, electron diffraction, light scattering, ab initio calculations, DFT theory,
	AtomicPhysics                                               Category = "physics.atom-ph"    // Atomic and molecular structure, spectra, collisions, and data. Atoms and molecules in external fields. Molecular dynamics and coherent and optical control. Cold atoms and molecules. Cold collisions. Optical lattices.
	BiologicalPhysics                                           Category = "physics.bio-ph"     // Molecular biophysics, cellular biophysics, neurological biophysics, membrane biophysics, single-molecule biophysics, ecological biophysics, quantum phenomena in biological systems (quantum biophysics), theoretical biophysics, molecular dynamics/modeling and simulation, game theory, biomechanics, bioinformatics, microorganisms, virology, evolution, biophysical methods
	ChemicalPhysics                                             Category = "physics.chem-ph"    // Experimental, computational, and theoretical physics of atoms, molecules, and clusters - Classical and quantum description of states, processes, and dynamics; spectroscopy, electronic structure, conformations, reactions, interactions, and phases. Chemical thermodynamics. Disperse systems. High pressure chemistry. Solid state chemistry. Surface and interface chemistry.
	ClassicalPhysics                                            Category = "physics.class-ph"   // Newtonian and relativistic dynamics; many particle systems; planetary motions; chaos in classical dynamics. Maxwell's equations and dynamics of charged systems and electromagnetic forces in materials. Vibrating systems such as membranes and cantilevers; optomechanics.
	ComputationalPhysics                                        Category = "physics.comp-ph"    // All aspects of computational science applied to physics.
	DataAnalysisStatisticsAndProbability                        Category = "physics.data-an"    // Methods, software and hardware for physics data analysis: data processing and storage; measurement methodology; statistical and mathematical aspects such as parametrization and uncertainties.
	PhysicsEducation                                            Category = "physics.ed-ph"      // Report of results of a research study, laboratory experience, assessment or classroom practice that represents a way to improve teaching and learning in physics. Also, report on misconceptions of students, textbook errors, and other similar information relative to promoting physics understanding.
	FluidDynamics                                               Category = "physics.flu-dyn"    // Turbulence, instabilities, incompressible/compressible flows, reacting flows. Aero/hydrodynamics, fluid-structure interactions, acoustics. Biological fluid dynamics, micro/nanofluidics, interfacial phenomena. Complex fluids, suspensions and granular flows, porous media flows.
	GeneralPhysics                                              Category = "physics.gen-ph"     // General Physics
	Geophysics                                                  Category = "physics.geo-ph"     // Atmospheric physics. Biogeosciences. Computational geophysics. Geographic location. Geoinformatics. Geophysical techniques. Hydrospheric geophysics. Magnetospheric physics. Mathematical geophysics. Planetology. Solar system. Solid earth geophysics. Space plasma physics. Mineral physics. High pressure physics.
	HistoryOfPhysics                                            Category = "physics.hist-ph"    // History and philosophy of all branches of physics, astrophysics, and cosmology, including appreciations of physicists.
	InstrumentationAndDetectors                                 Category = "physics.ins-det"    // Instrumentation and Detectors for research in natural science, including optical, molecular, atomic, nuclear and particle physics instrumentation and the associated electronics, services, infrastructure and control equipment.
	MedicalPhysics                                              Category = "physics.med-ph"     // Radiation therapy. Radiation dosimetry. Biomedical imaging modelling. Reconstruction, processing, and analysis. Biomedical system modelling and analysis. Health physics. New imaging or therapy modalities.
	Optics                                                      Category = "physics.optics"     // Adaptive optics. Astronomical optics. Atmospheric optics. Biomedical optics. Cardinal points. Collimation. Doppler effect. Fiber optics. Fourier optics. Geometrical optics (Gradient index optics. Holography. Infrared optics. Integrated optics. Laser applications. Laser optical systems
	PlasmaPhysics                                               Category = "physiscs.plasm-ph"  // Fundamental plasma physics. Magnetically Confined Plasmas (includes magnetic fusion energy research). High Energy Density Plasmas (inertial confinement plasmas, laser-plasma interactions).
	PopularPhysics                                              Category = "physics.pop-ph"     // Popular Physics
	PhysicsAndSociety                                           Category = "physics.soc-ph"     // Structure, dynamics and collective behavior of societies and groups (human or otherwise). Quantitative analysis of social networks and other complex networks. Physics and engineering of infrastructure and systems of broad societal impact (e.g., energy grids, transportation networks).
	SpacePhysics                                                Category = "physics.space-ph"   // Space plasma physics. Heliophysics. Space weather. Planetary magnetospheres, ionospheres and magnetotail. Auroras. Interplanetary space. Cosmic rays. Synchrotron radiation. Radio astronomy.
	QuantumPhysics                                              Category = "quant-ph"
) // Physics

const ( // Quantitative Biology
	QuantitativeBiologyBiomolecules            Category = "q-bio.BM" // DNA, RNA, proteins, lipids, etc.; molecular structures and folding kinetics; molecular interactions; single-molecule manipulation.
	QuantitativeBiologyCellBehavior            Category = "q-bio.CB" // Cell-cell signaling and interaction; morphogenesis and development; apoptosis; bacterial conjugation; viral-host interaction; immunology
	QuantitativeBiologyGenomics                Category = "q-bio.GN" // DNA sequencing and assembly; gene and motif finding; RNA editing and alternative splicing; genomic structure and processes (replication, transcription, methylation, etc); mutational processes.
	QuantitativeBiologyMolecularNetworks       Category = "q-bio.MN" // Gene regulation, signal transduction, proteomics, metabolomics, gene and enzymatic networks
	QuantitativeBiologyNeuronsAndCognition     Category = "q-bio.NC" // Synapse, cortex, neuronal dynamics, neural network, sensorimotor control, behavior, attention
	QuantitativeBiologyOther                   Category = "q-bio.OT" // Work in quantitative biology that does not fit into the other q-bio classifications
	QuantitativeBiologyPopulationsAndEvolution Category = "q-bio.PE" // Population dynamics, spatio-temporal and epidemiological models, dynamic speciation, co-evolution, biodiversity, foodwebs, aging; molecular evolution and phylogeny; directed evolution; origin of life
	QuantitativeBiologyQuantitativeMethods     Category = "q-bio.QM" // All experimental, numerical, statistical and mathematical contributions of value to biology
	QuantitativeBiologySubcellularProcesses    Category = "q-bio.SC" // Assembly and control of subcellular structures (channels, organelles, cytoskeletons, capsules, etc.); molecular motors, transport, subcellular localization; mitosis and meiosis
	QuantitativeBiologyTissuesAndOrgans        Category = "q-bio.TO" // Blood flow in vessels, biomechanics of bones, electrical waves, endocrine system, tumor growth
) // Quantitative Biology

const ( // Quantitative Finance
	ComputationalFinance           Category = "q-fin.CP" // Computational methods, including Monte Carlo, PDE, lattice and other numerical methods with applications to financial modeling
	Economics                      Category = "q-fin.EC" // q-fin.EC is an alias for econ.GN. Economics, including micro and macro economics, international economics, theory of the firm, labor economics, and other economic topics outside finance
	GeneralFinance                 Category = "q-fin.GN" // Development of general quantitative methodologies with applications in finance
	MathematicalFinance            Category = "q-fin.MF" // Mathematical and analytical methods of finance, including stochastic, probabilistic and functional analysis, algebraic, geometric and other methods
	PortfolioManagement            Category = "q-fin.PM" // Security selection and optimization, capital allocation, investment strategies and performance measurement
	PricingOfSecurities            Category = "q-fin.PR" // Valuation and hedging of financial securities, their derivatives, and structured products
	RiskManagement                 Category = "q-fin.RM" // Measurement and management of financial risks in trading, banking, insurance, corporate and other applications
	StatisticalFinance             Category = "q-fin.ST" // Statistical, econometric and econophysics analyses with applications to financial markets and economic data
	TradingAndMarketMicrostructure Category = "q-fin.TR" // Market microstructure, liquidity, exchange and auction design, automated trading, agent-based modeling and market-making
) // Quantitative Finance

const ( // Statistics
	StatisticsApplications    Category = "stat.AP" // Biology, Education, Epidemiology, Engineering, Environmental Sciences, Medical, Physical Sciences, Quality Control, Social Sciences
	StatisticsComputation     Category = "stat.CO" // Algorithms, Simulation, Visualization
	StatisticsMethodology     Category = "stat.ME" // Design, Surveys, Model Selection, Multiple Testing, Multivariate Methods, Signal and Image Processing, Time Series, Smoothing, Spatial Statistics, Survival Analysis, Nonparametric and Semiparametric Methods
	StatisticsMachineLearning Category = "stat.ML" // Covers machine learning papers (supervised, unsupervised, semi-supervised learning, graphical models, reinforcement learning, bandits, high dimensional inference, etc.) with a statistical or theoretical grounding
	StatisticsOther           Category = "stat.OT" // Work in statistics that does not fit into the other stat classifications
	StatisticsTheory          Category = "stat.TH" // stat.TH is an alias for math.ST. Asymptotics, Bayesian Inference, Decision Theory, Estimation, Foundations, Inference, Testing.
) // Statistics

func getLookupTable() map[Category]string {
	categoryMap := map[Category]string{
		HardwareArchitecture:                      "HardwareArchitecture",
		ArtificialIntelligence:                    "ArtificialIntelligence",
		ComputationAndLanguage:                    "ComputationAndLanguage",
		ComputationalComplexity:                   "ComputationalComplexity",
		ComputationalEngineeringFinanceAndScience: "ComputationalEngineeringFinanceAndScience",
		ComputationalGeometry:                     "ComputationalGeometry",
		ComputerScienceAndGameTheory:              "ComputerScienceAndGameTheory",
		ComputerVisionAndPatternRecognition:       "ComputerVisionAndPatternRecognition",
		ComputersAndSociety:                       "ComputersAndSociety",
		CryptographyAndSecurity:                   "CryptographyAndSecurity",
		Databases:                                 "Databases",
		DigitalLibraries:                          "DigitalLibraries",
		DiscreteMathematics:                       "DiscreteMathematics",
		DistributedParallelAndClusterComputing:    "DistributedParallelAndClusterComputing",
		EmergingTechnologies:                      "EmergingTechnologies",
		GeneralLiterature:                         "GeneralLiterature",
		Graphics:                                  "Graphics",
		HumanComputerInteraction:                  "HumanComputerInteraction",
		InformationRetrieval:                      "InformationRetrieval",
		InformationTheory:                         "InformationTheory",
		MachineLearning:                           "MachineLearning",
		LogicInComputerScience:                    "LogicInComputerScience",
		MathematicalSoftware:                      "MathematicalSoftware",
		MultiagentSystems:                         "MultiagentSystems",
		Multimedia:                                "Multimedia",
		NetworkAndInternetArchitecture:            "NetworkAndInternetArchitecture",
		NeuralAndEvolutionaryComputing:            "NeuralAndEvolutionaryComputing",
		NumericalAnalysis:                         "NumericalAnalysis",
		OperatingSystems:                          "OperatingSystems",
		CSOther:                                   "CSOther",
		Performance:                               "Performance",
		ProgrammingLanguages:                      "ProgrammingLanguages",
		Robotics:                                  "Robotics",
		SoftwareEngineering:                       "SoftwareEngineering",
		Sound:                                     "Sound",
		SymbolicComputation:                       "SymbolicComputation",
		SocialAndInformationNetworks:              "SocialAndInformationNetworks",
		CSSystemsAndControl:                       "SystemsAndControl", // Alias for eess.SY / SystemsAndControl

		// Economics (econ)
		Econometrics:         "Econometrics",
		GeneralEconomics:     "GeneralEconomics",
		TheoreticalEconomics: "TheoreticalEconomics",

		// Electrical Engineering and Systems Science (eess)
		AudioSndSpeechProcessing: "AudioSndSpeechProcessing",
		ImageAndVideoProcessing:  "ImageAndVideoProcessing",
		SignalProcessing:         "SignalProcessing",
		SystemsAndControl:        "SystemsAndControl",

		// Mathematics
		AlgebraicGeometry:                "AlgebraicGeometry",
		AlgebraicTopology:                "AlgebraicTopology",
		AnalysisOfPDEs:                   "AnalysisOfPDEs",
		CategoryTheory:                   "CategoryTheory",
		ClassicalAnalysisAndODEs:         "ClassicalAnalysisAndODEs",
		Combinatorics:                    "Combinatorics",
		CommutativeAlgebra:               "CommutativeAlgebra",
		ComplexVariables:                 "ComplexVariables",
		DataStructuresAndAlgorithms:      "DataStructuresAndAlgorithms",
		DifferentialGeometry:             "DifferentialGeometry",
		DynamicalSystems:                 "DynamicalSystems",
		FunctionalAnalysis:               "FunctionalAnalysis",
		FormalLanguagesAndAutomataTheory: "FormalLanguagesAndAutomataTheory",
		GeneralMathematics:               "GeneralMathematics",
		GeneralTopology:                  "GeneralTopology",
		GeometricTopology:                "GeometricTopology",
		GroupTheory:                      "GroupTheory",
		MathsHistoryAndOverview:          "MathsHistoryAndOverview",
		MathsInformationTheory:           "MathsInformationTheory",
		KTheoryAndHomology:               "KTheoryAndHomology",
		MathsLogic:                       "MathsLogic",
		MathsMathematicalPhysics:         "MathsMathematicalPhysics",
		MetricGeometry:                   "MetricGeometry",
		NumberTheory:                     "NumberTheory",
		MathsNumericalAnalysis:           "MathsNumericalAnalysis",
		OperatorAlgebras:                 "OperatorAlgebras",
		MathsOptimizationAndControl:      "MathsOptimizationAndControl",
		Probability:                      "Probability",
		QuantumAlgebra:                   "QuantumAlgebra",
		RepresentationTheory:             "RepresentationTheory",
		RingsAndAlgebra:                  "RingsAndAlgebra",
		MathsSpectralTheory:              "MathsSpectralTheory",
		SymplecticGeometry:               "SymplecticGeometry",
		MathsStatics:                     "MathsStatics",

		// Physics
		Physics:                                  "Physics",
		Astrophysics:                             "Astrophysics",
		CosmologyAndNongalacticAstrophysics:      "CosmologyAndNongalacticAstrophysics",
		EarthAndPlanetaryAstrophysics:            "EarthAndPlanetaryAstrophysics",
		AstrophysicsOfGalaxies:                   "AstrophysicsOfGalaxies",
		HighEnergyAstrophysicalPhenomena:         "HighEnergyAstrophysicalPhenomena",
		InstrumentationAndMethodsForAstrophysics: "InstrumentationAndMethodsForAstrophysics",
		SolarAndStellarAstrophysics:              "SolarAndStellarAstrophysics",
		CondensedMatter:                          "CondensedMatter",
		DisorderedSystemsAndNeuralNetworks:       "DisorderedSystemsAndNeuralNetworks",
		MesoscaleAndNanoscalePhysics:             "MesoscaleAndNanoscalePhysics",
		MaterialsScience:                         "MaterialsScience",
		PhysicsOther:                             "PhysicsOther",
		QuantumGases:                             "QuantumGases",
		SoftCondensedMatter:                      "SoftCondensedMatter",
		StatisticalMechanics:                     "StatisticalMechanics",
		StronglyCorrelatedElectrons:              "StronglyCorrelatedElectrons",
		Superconductivity:                        "Superconductivity",
		GeneralRelativityAndQuantumCosmology:     "GeneralRelativityAndQuantumCosmology",
		HighEnergyPhysic:                         "HighEnergyPhysic",
		HighEnergyPhysicsExperiment:              "HighEnergyPhysicsExperiment",
		HighEnergyPhysicsLattice:                 "HighEnergyPhysicsLattice",
		HighEnergyPhysicsPhenomenology:           "HighEnergyPhysicsPhenomenology",
		HighEnergyPhysicsTheory:                  "HighEnergyPhysicsTheory",
		MathematicalPhysics:                      "MathematicalPhysics",
		NonLinearSciences:                        "NonLinearSciences",
		NonLinearSciencesAdaptationAndSelfOrganizingSystemsCategory: "NonLinearSciencesAdaptationAndSelfOrganizingSystemsCategory",
		NonLinearSciencesCellularAutomataAndLattice:                 "NonLinearSciencesCellularAutomataAndLattice",
		NonLinearSciencesChaoticDynamics:                            "NonLinearSciencesChaoticDynamics",
		ExactlySolvableAndIntegrableSystems:                         "ExactlySolvableAndIntegrableSystems",
		PatternFormationAndSolutions:                                "PatternFormationAndSolutions",
		NuclearExperiment:                                           "NuclearExperiment",
		NuclearTheory:                                               "NuclearTheory",
		AcceleratorPhysics:                                          "AcceleratorPhysics",
		AtmoshpericAndOceanicPhysics:                                "AtmoshpericAndOceanicPhysics",
		AtomicPhysics:                                               "AtomicPhysics",
		AppliedPhysics:                                              "AppliedPhysics",
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

		// Quantitative Biology
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

		// Quantitative Finance
		ComputationalFinance:           "Computational Finance",
		Economics:                      "Economics",
		GeneralFinance:                 "General Finance",
		MathematicalFinance:            "Mathematical Finance",
		PortfolioManagement:            "Portfolio Management",
		PricingOfSecurities:            "Pricing of Securities",
		RiskManagement:                 "Risk Management",
		StatisticalFinance:             "Statistical Finance",
		TradingAndMarketMicrostructure: "Trading and Market Microstructure",

		// Statistics
		StatisticsApplications:    "StatisticsApplications",
		StatisticsComputation:     "StatisticsComputation",
		StatisticsMachineLearning: "StatisticsMachineLearning",
		StatisticsMethodology:     "StatisticsMethodology",
		StatisticsTheory:          "StatisticsTheory",
		StatisticsOther:           "StatisticsOther",
	}
	return categoryMap
}
