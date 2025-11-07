use std::num::NonZeroU8;

// Example IDs:
// NCIT:C25330
// NCIT:R100
// NCIT:P383
// MS:1000014
// UO:0000245
// BAO_0000925
// BFO:0000015
// GAZ:00002933
// AfPO_0000233 // In URL
// BTO:0004947
// PRIDE:0000521
// MOD:01188
// XLMOD:07097
// GNO:G00001NT // Name and ID are identical
// GNO:00000015

/// Controlled vocabularies all ontobee listed controlled vocabulaires are present as well as PRIDE and RESID
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(u16)]
#[allow(
    clippy::upper_case_acronyms,
    non_camel_case_types,
    clippy::doc_markdown
)]
pub enum ControlledVocabulary {
    ///  Alzheimer's Disease Ontology
    ADO,
    ///  Anatomical Entity Ontology
    AEO,
    ///  Allotrope Foundation Ontology
    AFO,
    ///  African Population Ontology
    AfPO,
    ///  Agronomy Ontology
    AGRO,
    ///  Ontology for the Anatomy of the Insect SkeletoMuscular system (AISM)
    AISM,
    ///  The Amphioxus Development and Anatomy Ontology
    AMPHX,
    ///  Ascomycete phenotype ontology
    APO,
    ///  Apollo Structured Vocabulary
    APOLLO_SV,
    ///  Antibiotic Resistance Ontology
    ARO,
    ///  BioAssay Ontology
    BAO,
    ///  Beta Cell Genomics Ontology
    BCGO,
    ///  Biological Collections Ontology
    BCO,
    ///  Basic Formal Ontology
    BFO,
    ///  Basic Formal Ontology (BFO) 1.1
    BFO11,
    ///  Bacteria Infectious Disease Ontology
    BIDO,
    ///  Biomarker Ontology
    BMONT,
    ///  Biological Spatial Ontology
    BSPO,
    ///  BRENDA tissue / enzyme source
    BTO,
    ///  Common Anatomy Reference Ontology
    CARO,
    ///  Common Coordinate Framework Ontology
    CCFO,
    ///  Comparative Data Analysis Ontology
    CDAO,
    ///  Compositional Dietary Nutrition Ontology
    CDNO,
    ///  Cephalopod Ontology
    CEPH,
    ///  Chemical Entities of Biological Interest
    CHEBI,
    ///  Chemical Information Ontology
    CHEMINF,
    ///  CHEBI Integrated Role Ontology
    CHIRO,
    ///  Chemical Methods Ontology
    CHMO,
    ///  Coronavirus Infectious Disease Ontology
    CIDO,
    ///  Confidence Information Ontology
    CIO,
    ///  Cell Ontology
    CL,
    ///  Collembola Anatomy Ontology
    CLAO,
    ///  Cell Line Ontology
    CLO,
    ///  Cell Line Ontology - Human Pluripotent Stem Cell Registry
    CLO_hPSCreg,
    ///  Cell Line Ontology - Chinese National Infrastructure of Cell Line Resource
    CLO_NICR,
    ///  Clytia hemisphaerica Development and Anatomy Ontology
    CLYH,
    ///  CranioMaxilloFacial ontology
    CMF,
    ///  Clinical measurement ontology
    CMO,
    ///  Core Ontology for Biology and Biomedicine
    COB,
    ///  Coleoptera Anatomy Ontology (COLAO)
    COLAO,
    ///  Contributor Role Ontology
    CRO,
    ///  OAE CTCAE view
    CTCAE_OAEview,
    ///  Ctenophore Ontology
    CTENO,
    ///  CTO: Core Ontology of Clinical Trials
    CTO,
    ///  Cardiovascular Disease Ontology
    CVDO,
    ///  Ontology of Document Acts
    d_acts,
    ///  Dictyostelium discoideum anatomy
    DDANAT,
    ///  Dictyostelium discoideum phenotype ontology
    DDPHENO,
    ///  Drug-drug Interaction and Drug-drug Interaction Evidence Ontology
    DIDEO,
    ///  The Drug-Drug Interactions Ontology
    DINTO,
    ///  Disease Drivers Ontology
    DISDRIV,
    ///  Human Disease Ontology
    DOID,
    ///  The Drug Ontology
    DRON,
    ///  Data Use Ontology
    DUO,
    ///  The Echinoderm Anatomy and Development Ontology
    ECAO,
    ///  Evidence and Conclusion Ontology
    ECO,
    ///  An ontology of core ecological entities
    ECOCORE,
    ///  Environmental conditions, treatments and exposures ontology
    ECTO,
    ///  EMBRACE Data and Methods
    EDAM,
    ///  Experimental Factor Ontology
    EFO,
    ///  Human developmental anatomy, abstract
    EHDAA2,
    ///  Mouse Developmental Anatomy Ontology
    EMAPA,
    ///  Environment Ontology
    ENVO,
    ///  Epilepsy Ontology
    EPIO,
    ///  Epidemiology Ontology
    EPO,
    ///  eagle-i resource ontology
    ERO,
    ///  VEuPathDB ontology
    EUPATH,
    ///  Exercise Medicine Ontology
    EXMO,
    ///  Exposure ontology
    ExO,
    ///  Fungal gross anatomy
    FAO,
    ///  Biological Imaging Methods Ontology
    FBbi,
    ///  Drosophila gross anatomy
    FBbt,
    ///  Drosophila Phenotype Ontology
    FBcv,
    ///  Drosophila development
    FBdv,
    ///  Food Interactions with Drugs Evidence Ontology
    FIDEO,
    ///  Physico-chemical methods and properties
    FIX,
    ///  Flora Phenotype Ontology
    FLOPO,
    ///  Foundational Model of Anatomy Ontology (subset)
    FMA,
    ///  Food-Biomarker Ontology
    FOBI,
    ///  Food Ontology
    FOODON,
    ///  FuTRES Ontology of Vertebrate Traits
    FOVT,
    ///  Fission Yeast Phenotype Ontology
    FYPO,
    ///  Plant Gall Ontology
    GALLONT,
    ///  Gazetteer
    GAZ,
    ///  Genomics Cohorts Knowledge Ontology
    GECKO,
    ///  Genomic Epidemiology Ontology
    GENEPIO,
    ///  Genotype Ontology
    GENO,
    ///  Geographical Entity Ontology
    GEO,
    ///  Glycan Naming and Subsumption Ontology (GNOme)
    GNO,
    ///  Gene Ontology
    GO,
    ///  Gender, Sex, and Sexual Orientation (GSSO) ontology
    GSSO,
    ///  Human Ancestry Ontology
    HANCESTRO,
    ///  Hymenoptera Anatomy Ontology
    HAO,
    ///  Human Interaction Network Ontology
    HINO,
    ///  Homology Ontology
    HOM,
    ///  Human Phenotype Ontology (HPO)
    HP,
    ///  Human Developmental Stages
    HsapDv,
    ///  Health Surveillance Ontology
    HSO,
    ///  Hypertension Ontology
    HTN,
    ///  Information Artifact Ontology
    IAO,
    ///  IAO Ontology Metadata
    IAO_Onto_Meta,
    ///  International Classification of Disease Ontology
    ICDO,
    ///  Integrative and Conjugative Element Ontology
    ICEO,
    ///  Informed Consent Ontology
    ICO,
    ///  Infectious Disease Ontology
    IDO,
    ///  Brucellosis Ontology
    IDOBRU,
    ///  Malaria Ontology
    IDOMAL,
    ///  Interaction Network Ontology
    INO,
    ///  Kinetic Simulation Algorithm Ontology
    KISAO,
    ///  Kidney Tissue Atlas Ontology
    KTAO,
    ///  clinical LABoratory Ontology
    LABO,
    ///  Lepidoptera Anatomy Ontology
    LEPAO,
    ///  CLO LINCS view
    LINCS_CLOview,
    ///  LTHIDO
    LTHIDO,
    ///  Mouse adult gross anatomy
    MA,
    ///  Mathematical modeling ontology
    MAMO,
    ///  Medical Action Ontology
    MAXO,
    ///  Microbial Conditions Ontology
    MCO,
    ///  Model Card Report Ontology
    MCRO,
    ///  Mental Functioning Ontology
    MF,
    ///  Mammalian Feeding Muscle Ontology
    MFMO,
    ///  Emotion Ontology
    MFOEM,
    ///  Mental Disease Ontology
    MFOMD,
    ///  Molecular Interactions Controlled Vocabulary
    MI,
    ///  MIAPA Ontology
    MIAPA,
    ///  Ontology of Prokaryotic Phenotypic and Metabolic Characters
    MICRO,
    ///  Mycosis Infectious Disease Ontology
    MIDO,
    ///  microRNA Ontology
    miRNAO,
    ///  Mosquito insecticide resistance
    MIRO,
    ///  Measurement method ontology
    MMO,
    ///  Mouse Developmental Stages
    MmusDv,
    ///  Protein modification
    MOD,
    ///  Mondo Disease Ontology
    MONDO,
    ///  Molecular Process Ontology
    MOP,
    ///  Mammalian Phenotype Ontology
    MP,
    ///  Mouse pathology
    MPATH,
    ///  Minimum PDDI Information Ontology
    MPIO,
    ///  MHC Restriction Ontology
    MRO,
    ///  The PSI-MS Controlled Vocabulary [https://www.ebi.ac.uk/ols4/ontologies/ms](https://www.ebi.ac.uk/ols4/ontologies/ms)
    MS,
    ///  Neuro Behavior Ontology
    NBO,
    ///  NCBI organismal classification
    NCBITaxon,
    ///  NCI Thesaurus OBO Edition
    NCIT,
    ///  Non-Coding RNA Ontology
    NCRO,
    ///  National Drug File Reference Terminology
    NDF_RT,
    ///  NFO
    NFO,
    ///  Next Generation Biobanking Ontology
    NGBO,
    ///  NOMEN - A nomenclatural ontology for biological names
    NOMEN,
    ///  Ontology of Adverse Events
    OAE,
    ///  Ontology of Arthropod Circulatory Systems
    OARCS,
    ///  Ontology of Biological Attributes
    OBA,
    ///  Ontology of Biological and Clinical Statistics
    OBCS,
    ///  Ontology for Biomedical Investigations
    OBI,
    ///  OBI NIAID-GSC-BRC view
    OBI_NIAID_GSC_BRC_view,
    ///  Ontology for Biobanking
    OBIB,
    ///  OBI web service, development version
    OBIws,
    ///  Occupation Ontology
    OCCO,
    ///  Ontology of Chemical Elements
    OCE,
    ///  Ontology of Chinese Medicine for Rheumatism
    OCMR,
    ///  OCRe
    OCRe,
    ///  Ontology of Cardiovascular Drug Adverse Events
    OCVDAE,
    ///  Ontology of Drug Adverse Events
    ODAE,
    ///  Ontologyof Drug Neuropathy Adverse Events
    ODNAE,
    ///  The Ontology of Genes and Genomes
    OGG,
    ///  Ontology of Genes and Genomes - Arabidopsis thaliana
    OGG_At,
    ///  Ontology of Genes and Genomes - Brucella
    OGG_Bru,
    ///  Ontology of Genes and Genomes - Caenorhabditis elegans
    OGG_Ce,
    ///  Ontology of Genes and Genomes - Fruit Fly
    OGG_Dm,
    ///  Ontology of Genes and Genomes - Zebrafish
    OGG_Dr,
    ///  Ontology of Genes and Genomes - Mouse
    OGG_Mm,
    ///  Ontology of Genes and Genomes - Plasmodium falciparum
    OGG_Pf,
    ///  Ontology of Genes and Genomes - Yeast
    OGG_Sc,
    ///  Ontology for genetic interval
    OGI,
    ///  Ontology for General Medical Science
    OGMS,
    ///  Ontology of Genetic Susceptibility Factor
    OGSF,
    ///  Oral Health and Disease Ontology
    OHD,
    ///  Ontology of Host-Microbiome Interactions
    OHMI,
    ///  Ontology of Host Pathogen Interactions
    OHPI,
    ///  Medaka Developmental Stages
    OlatDv,
    ///  Ontology for Linking Biological and Medical Ontologies
    OLOBO,
    ///  Ontology of units of Measure
    OM,
    ///  Ontologized MIABIS
    OMIABIS,
    ///  Ontology for MIRNA Target
    OMIT,
    ///  OBO Metadata Ontology
    OMO,
    ///  Ontology of Microbial Phenotypes
    OMP,
    ///  Ontology for Modeling and Representation of Social Entities
    OMRSE,
    ///  Ontology for Nutritional Epidemiology
    ONE,
    ///  Ontology for Nutritional Studies
    ONS,
    ///  OntoAvida: ontology for Avida digital evolution platform
    ONTOAVIDA,
    ///  Obstetric and Neonatal Ontology
    ONTONEO,
    ///  Ontology of Organizational Structures of Trauma centers and Trauma systems
    OOSTT,
    ///  Ontology for Parasite LifeCycle
    OPL,
    ///  Ontology of Precision Medicine and Investigation
    OPMI,
    ///  Ontology of RNA Sequencing
    ORNASEQ,
    ///  Ontology for Stem Cell Investigations
    OSCI,
    ///  Ontology of Vaccine Adverse Events
    OVAE,
    ///  Phenotype And Trait Ontology
    PATO,
    ///  Provisional Cell Ontology
    PCL,
    ///  Population and Community Ontology
    PCO,
    ///  The Prescription of Drugs Ontology
    PDRO,
    ///  Platynereis Developmental Stages
    PdumDv,
    ///  Plant Experimental Conditions Ontology
    PECO,
    ///  Promoting Health Aging through Semantic Enrichment of Solitude Research (PHASES) Ontology
    PHASES,
    ///  Pathogen Host Interaction Phenotype Ontology
    PHIPO,
    ///  Parasite Infectious Disease Ontology
    PIDO,
    ///  planaria-ontology
    PLANA,
    ///  Planarian Phenotype Ontology
    PLANP,
    ///  Proper Name Ontology
    PNO,
    ///  Plant Ontology
    PO,
    ///  Porifera Ontology
    PORO,
    ///  Plant Phenology Ontology
    PPO,
    ///  PRotein Ontology (PRO)
    PR,
    /// The PRIDE Controlled Vocabulary <https://www.ebi.ac.uk/ols4/ontologies/pride>
    PRIDE,
    ///  Process Chemistry Ontology
    PROCO,
    ///  Performance Summary Display Ontology
    PSDO,
    ///  Plant Stress Ontology
    PSO,
    ///  Pathway ontology
    PW,
    ///  Radiation Biology Ontology
    RBO,
    /// Reagent Ontology
    REO,
    /// RESID top down modifications ontology
    RESID,
    ///  Physico-chemical process
    REX,
    ///  RNA ontology
    RNAO,
    ///  Relation Ontology
    RO,
    ///  Rat Strain Ontology
    RS,
    ///  Name Reaction Ontology
    RXNO,
    ///  Systems Biology Ontology
    SBO,
    ///  Sickle Cell Disease Ontology
    SCDO,
    ///  Sustainability Core Ontology
    SCO,
    ///  Sustainable Development Goals Interface Ontology
    SDGIO,
    ///  Sample processing and separation techniques
    SEP,
    ///  Scientific Evidence and Provenance Information Ontology
    SEPIO,
    ///  Social Insect Behavior Ontology
    SIBO,
    ///  Semanticscience Integrated Ontology
    SIO,
    ///  Space Life Sciences Ontology
    SLSO,
    ///  Sequence types and features ontology
    SO,
    ///  Spider Ontology
    SPD,
    ///  The Statistical Methods Ontology
    STATO,
    ///  Software ontology
    SWO,
    ///  Symptom Ontology
    SYMP,
    ///  terms4FAIRskills
    T4FS,
    ///  Tick Anatomy Ontology
    TADS,
    ///  Taxonomic rank vocabulary
    TAXRANK,
    ///  Mosquito gross anatomy ontology
    TGMA,
    ///  Plant Trait Ontology
    TO,
    ///  Pathogen Transmission Ontology
    TRANS,
    ///  Transportation System Ontology
    TSO,
    ///  Teleost taxonomy ontology
    TTO,
    ///  Toxic Process Ontology
    TXPO,
    ///  Uberon multi-species anatomy ontology
    UBERON,
    ///  Uber anatomy ontology, basic version
    uberon_basic,
    ///  uberon_collected_metazoa
    uberon_collected_metazoa,
    ///  uberon_depictions
    uberon_depictions,
    ///  Units of measurement ontology
    UO,
    ///  Unipathway
    UPA,
    ///  Unified Phenotype Ontology (uPheno)
    UPHENO,
    ///  Variation Ontology
    VariO,
    ///  Vertebrate Breed Ontology
    VBO,
    ///  Vaccination Informed Consent Ontology
    VICO,
    ///  Virus Infectious Disease Ontology
    VIDO,
    ///  Vaccine Investigation Ontology
    VIO,
    ///  VIVO-ISF
    VIVO_ISF,
    ///  Vaccine Ontology
    VO,
    ///  Vertebrate trait ontology
    VT,
    ///  Vertebrate Taxonomy Ontology
    VTO,
    ///  C. elegans Gross Anatomy Ontology
    WBbt,
    ///  C. elegans development ontology
    WBls,
    ///  C. elegans phenotype
    WBPhenotype,
    ///  Xenopus Anatomy Ontology
    XAO,
    ///  Experimental condition ontology
    XCO,
    ///  Cross-linker reagents ontology
    XL,
    ///  HUPO-PSI cross-linking and derivatization reagents controlled vocabulary
    XLMOD,
    ///  Xenopus Phenotype Ontology
    XPO,
    ///  Zebrafish Experimental Conditions Ontology
    ZECO,
    ///  Zebrafish anatomy and development ontology
    ZFA,
    ///  Zebrafish developmental stages ontology
    ZFS,
    ///  Zebrafish Phenotype Ontology
    ZP,
    /// An unknown ontology
    Unknown,
}

impl std::str::FromStr for ControlledVocabulary {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ADO" => Ok(Self::ADO),
            "AEO" => Ok(Self::AEO),
            "AFO" => Ok(Self::AFO),
            "AfPO" => Ok(Self::AfPO),
            "AGRO" => Ok(Self::AGRO),
            "AISM" => Ok(Self::AISM),
            "AMPHX" => Ok(Self::AMPHX),
            "APO" => Ok(Self::APO),
            "APOLLO_SV" => Ok(Self::APOLLO_SV),
            "ARO" => Ok(Self::ARO),
            "BAO" => Ok(Self::BAO),
            "BCGO" => Ok(Self::BCGO),
            "BCO" => Ok(Self::BCO),
            "BFO" => Ok(Self::BFO),
            "BFO11" => Ok(Self::BFO11),
            "BIDO" => Ok(Self::BIDO),
            "BMONT" => Ok(Self::BMONT),
            "BSPO" => Ok(Self::BSPO),
            "BTO" => Ok(Self::BTO),
            "CARO" => Ok(Self::CARO),
            "CCFO" => Ok(Self::CCFO),
            "CDAO" => Ok(Self::CDAO),
            "CDNO" => Ok(Self::CDNO),
            "CEPH" => Ok(Self::CEPH),
            "CHEBI" => Ok(Self::CHEBI),
            "CHEMINF" => Ok(Self::CHEMINF),
            "CHIRO" => Ok(Self::CHIRO),
            "CHMO" => Ok(Self::CHMO),
            "CIDO" => Ok(Self::CIDO),
            "CIO" => Ok(Self::CIO),
            "CL" => Ok(Self::CL),
            "CLAO" => Ok(Self::CLAO),
            "CLO" => Ok(Self::CLO),
            "CLO-hPSCreg" | "CLO_hPSCreg" => Ok(Self::CLO_hPSCreg),
            "CLO-NICR" | "CLO_NICR" => Ok(Self::CLO_NICR),
            "CLYH" => Ok(Self::CLYH),
            "CMF" => Ok(Self::CMF),
            "CMO" => Ok(Self::CMO),
            "COB" => Ok(Self::COB),
            "COLAO" => Ok(Self::COLAO),
            "CRO" => Ok(Self::CRO),
            "CTCAE-OAEview" | "CTCAE_OAEview" => Ok(Self::CTCAE_OAEview),
            "CTENO" => Ok(Self::CTENO),
            "CTO" => Ok(Self::CTO),
            "CVDO" => Ok(Self::CVDO),
            "d-acts" | "d_acts" => Ok(Self::d_acts),
            "DDANAT" => Ok(Self::DDANAT),
            "DDPHENO" => Ok(Self::DDPHENO),
            "DIDEO" => Ok(Self::DIDEO),
            "DINTO" => Ok(Self::DINTO),
            "DISDRIV" => Ok(Self::DISDRIV),
            "DOID" => Ok(Self::DOID),
            "DRON" => Ok(Self::DRON),
            "DUO" => Ok(Self::DUO),
            "ECAO" => Ok(Self::ECAO),
            "ECO" => Ok(Self::ECO),
            "ECOCORE" => Ok(Self::ECOCORE),
            "ECTO" => Ok(Self::ECTO),
            "EDAM" => Ok(Self::EDAM),
            "EFO" => Ok(Self::EFO),
            "EHDAA2" => Ok(Self::EHDAA2),
            "EMAPA" => Ok(Self::EMAPA),
            "ENVO" => Ok(Self::ENVO),
            "EPIO" => Ok(Self::EPIO),
            "EPO" => Ok(Self::EPO),
            "ERO" => Ok(Self::ERO),
            "EUPATH" => Ok(Self::EUPATH),
            "EXMO" => Ok(Self::EXMO),
            "ExO" => Ok(Self::ExO),
            "FAO" => Ok(Self::FAO),
            "FBbi" => Ok(Self::FBbi),
            "FBbt" => Ok(Self::FBbt),
            "FBcv" => Ok(Self::FBcv),
            "FBdv" => Ok(Self::FBdv),
            "FIDEO" => Ok(Self::FIDEO),
            "FIX" => Ok(Self::FIX),
            "FLOPO" => Ok(Self::FLOPO),
            "FMA" => Ok(Self::FMA),
            "FOBI" => Ok(Self::FOBI),
            "FOODON" => Ok(Self::FOODON),
            "FOVT" => Ok(Self::FOVT),
            "FYPO" => Ok(Self::FYPO),
            "GALLONT" => Ok(Self::GALLONT),
            "GAZ" => Ok(Self::GAZ),
            "GECKO" => Ok(Self::GECKO),
            "GENEPIO" => Ok(Self::GENEPIO),
            "GENO" => Ok(Self::GENO),
            "GEO" => Ok(Self::GEO),
            "GNO" | "GNOME" | "G" => Ok(Self::GNO),
            "GO" => Ok(Self::GO),
            "GSSO" => Ok(Self::GSSO),
            "HANCESTRO" => Ok(Self::HANCESTRO),
            "HAO" => Ok(Self::HAO),
            "HINO" => Ok(Self::HINO),
            "HOM" => Ok(Self::HOM),
            "HP" => Ok(Self::HP),
            "HsapDv" => Ok(Self::HsapDv),
            "HSO" => Ok(Self::HSO),
            "HTN" => Ok(Self::HTN),
            "IAO" => Ok(Self::IAO),
            "IAO-Onto-Meta" | "IAO_Onto_Meta" => Ok(Self::IAO_Onto_Meta),
            "ICDO" => Ok(Self::ICDO),
            "ICEO" => Ok(Self::ICEO),
            "ICO" => Ok(Self::ICO),
            "IDO" => Ok(Self::IDO),
            "IDOBRU" => Ok(Self::IDOBRU),
            "IDOMAL" => Ok(Self::IDOMAL),
            "INO" => Ok(Self::INO),
            "KISAO" => Ok(Self::KISAO),
            "KTAO" => Ok(Self::KTAO),
            "LABO" => Ok(Self::LABO),
            "LEPAO" => Ok(Self::LEPAO),
            "LINCS-CLOview" | "LINCS_CLOview" => Ok(Self::LINCS_CLOview),
            "LTHIDO" => Ok(Self::LTHIDO),
            "MA" => Ok(Self::MA),
            "MAMO" => Ok(Self::MAMO),
            "MAXO" => Ok(Self::MAXO),
            "MCO" => Ok(Self::MCO),
            "MCRO" => Ok(Self::MCRO),
            "MF" => Ok(Self::MF),
            "MFMO" => Ok(Self::MFMO),
            "MFOEM" => Ok(Self::MFOEM),
            "MFOMD" => Ok(Self::MFOMD),
            "MI" => Ok(Self::MI),
            "MIAPA" => Ok(Self::MIAPA),
            "MICRO" => Ok(Self::MICRO),
            "MIDO" => Ok(Self::MIDO),
            "miRNAO" => Ok(Self::miRNAO),
            "MIRO" => Ok(Self::MIRO),
            "MMO" => Ok(Self::MMO),
            "MmusDv" => Ok(Self::MmusDv),
            "MOD" | "PSI-MOD" => Ok(Self::MOD),
            "MONDO" => Ok(Self::MONDO),
            "MOP" => Ok(Self::MOP),
            "MP" => Ok(Self::MP),
            "MPATH" => Ok(Self::MPATH),
            "MPIO" => Ok(Self::MPIO),
            "MRO" => Ok(Self::MRO),
            "MS" | "PSI-MS" => Ok(Self::MS),
            "NBO" => Ok(Self::NBO),
            "NCBITaxon" => Ok(Self::NCBITaxon),
            "NCIT" => Ok(Self::NCIT),
            "NCRO" => Ok(Self::NCRO),
            "NDF-RT" | "NDF_RT" => Ok(Self::NDF_RT),
            "NFO" => Ok(Self::NFO),
            "NGBO" => Ok(Self::NGBO),
            "NOMEN" => Ok(Self::NOMEN),
            "OAE" => Ok(Self::OAE),
            "OARCS" => Ok(Self::OARCS),
            "OBA" => Ok(Self::OBA),
            "OBCS" => Ok(Self::OBCS),
            "OBI" => Ok(Self::OBI),
            "OBI-NIAID-GSC-BRC-view" | "OBI_NIAID_GSC_BRC_view" => Ok(Self::OBI_NIAID_GSC_BRC_view),
            "OBIB" => Ok(Self::OBIB),
            "OBIws" => Ok(Self::OBIws),
            "OCCO" => Ok(Self::OCCO),
            "OCE" => Ok(Self::OCE),
            "OCMR" => Ok(Self::OCMR),
            "OCRe" => Ok(Self::OCRe),
            "OCVDAE" => Ok(Self::OCVDAE),
            "ODAE" => Ok(Self::ODAE),
            "ODNAE" => Ok(Self::ODNAE),
            "OGG" => Ok(Self::OGG),
            "OGG-At" | "OGG_At" => Ok(Self::OGG_At),
            "OGG-Bru" | "OGG_Bru" => Ok(Self::OGG_Bru),
            "OGG-Ce" | "OGG_Ce" => Ok(Self::OGG_Ce),
            "OGG-Dm" | "OGG_Dm" => Ok(Self::OGG_Dm),
            "OGG-Dr" | "OGG_Dr" => Ok(Self::OGG_Dr),
            "OGG-Mm" | "OGG_Mm" => Ok(Self::OGG_Mm),
            "OGG-Pf" | "OGG_Pf" => Ok(Self::OGG_Pf),
            "OGG-Sc" | "OGG_Sc" => Ok(Self::OGG_Sc),
            "OGI" => Ok(Self::OGI),
            "OGMS" => Ok(Self::OGMS),
            "OGSF" => Ok(Self::OGSF),
            "OHD" => Ok(Self::OHD),
            "OHMI" => Ok(Self::OHMI),
            "OHPI" => Ok(Self::OHPI),
            "OlatDv" => Ok(Self::OlatDv),
            "OLOBO" => Ok(Self::OLOBO),
            "OM" => Ok(Self::OM),
            "OMIABIS" => Ok(Self::OMIABIS),
            "OMIT" => Ok(Self::OMIT),
            "OMO" => Ok(Self::OMO),
            "OMP" => Ok(Self::OMP),
            "OMRSE" => Ok(Self::OMRSE),
            "ONE" => Ok(Self::ONE),
            "ONS" => Ok(Self::ONS),
            "ONTOAVIDA" => Ok(Self::ONTOAVIDA),
            "ONTONEO" => Ok(Self::ONTONEO),
            "OOSTT" => Ok(Self::OOSTT),
            "OPL" => Ok(Self::OPL),
            "OPMI" => Ok(Self::OPMI),
            "ORNASEQ" => Ok(Self::ORNASEQ),
            "OSCI" => Ok(Self::OSCI),
            "OVAE" => Ok(Self::OVAE),
            "PATO" => Ok(Self::PATO),
            "PCL" => Ok(Self::PCL),
            "PCO" => Ok(Self::PCO),
            "PDRO" => Ok(Self::PDRO),
            "PdumDv" => Ok(Self::PdumDv),
            "PECO" => Ok(Self::PECO),
            "PHASES" => Ok(Self::PHASES),
            "PHIPO" => Ok(Self::PHIPO),
            "PIDO" => Ok(Self::PIDO),
            "PLANA" => Ok(Self::PLANA),
            "PLANP" => Ok(Self::PLANP),
            "PNO" => Ok(Self::PNO),
            "PO" => Ok(Self::PO),
            "PORO" => Ok(Self::PORO),
            "PPO" => Ok(Self::PPO),
            "PR" => Ok(Self::PR),
            "PRIDE" => Ok(Self::PRIDE),
            "PROCO" => Ok(Self::PROCO),
            "PSDO" => Ok(Self::PSDO),
            "PSO" => Ok(Self::PSO),
            "PW" => Ok(Self::PW),
            "RBO" => Ok(Self::RBO),
            "REO" => Ok(Self::REO),
            "RESID" => Ok(Self::RESID),
            "REX" => Ok(Self::REX),
            "RNAO" => Ok(Self::RNAO),
            "RO" => Ok(Self::RO),
            "RS" => Ok(Self::RS),
            "RXNO" => Ok(Self::RXNO),
            "SBO" => Ok(Self::SBO),
            "SCDO" => Ok(Self::SCDO),
            "SCO" => Ok(Self::SCO),
            "SDGIO" => Ok(Self::SDGIO),
            "SEP" => Ok(Self::SEP),
            "SEPIO" => Ok(Self::SEPIO),
            "SIBO" => Ok(Self::SIBO),
            "SIO" => Ok(Self::SIO),
            "SLSO" => Ok(Self::SLSO),
            "SO" => Ok(Self::SO),
            "SPD" => Ok(Self::SPD),
            "STATO" => Ok(Self::STATO),
            "SWO" => Ok(Self::SWO),
            "SYMP" => Ok(Self::SYMP),
            "T4FS" => Ok(Self::T4FS),
            "TADS" => Ok(Self::TADS),
            "TAXRANK" => Ok(Self::TAXRANK),
            "TGMA" => Ok(Self::TGMA),
            "TO" => Ok(Self::TO),
            "TRANS" => Ok(Self::TRANS),
            "TSO" => Ok(Self::TSO),
            "TTO" => Ok(Self::TTO),
            "TXPO" => Ok(Self::TXPO),
            "UBERON" => Ok(Self::UBERON),
            "uberon-basic" | "uberon_basic" => Ok(Self::uberon_basic),
            "uberon_collected_metazoa" => Ok(Self::uberon_collected_metazoa),
            "uberon_depictions" => Ok(Self::uberon_depictions),
            "UO" => Ok(Self::UO),
            "UPA" => Ok(Self::UPA),
            "UPHENO" => Ok(Self::UPHENO),
            "VariO" => Ok(Self::VariO),
            "VBO" => Ok(Self::VBO),
            "VICO" => Ok(Self::VICO),
            "VIDO" => Ok(Self::VIDO),
            "VIO" => Ok(Self::VIO),
            "VIVO-ISF" | "VIVO_ISF" => Ok(Self::VIVO_ISF),
            "VO" => Ok(Self::VO),
            "VT" => Ok(Self::VT),
            "VTO" => Ok(Self::VTO),
            "WBbt" => Ok(Self::WBbt),
            "WBls" => Ok(Self::WBls),
            "WBPhenotype" => Ok(Self::WBPhenotype),
            "XAO" => Ok(Self::XAO),
            "XCO" => Ok(Self::XCO),
            "XL" => Ok(Self::XL),
            "XLMOD" => Ok(Self::XLMOD),
            "XPO" => Ok(Self::XPO),
            "ZECO" => Ok(Self::ZECO),
            "ZFA" => Ok(Self::ZFA),
            "ZFS" => Ok(Self::ZFS),
            "ZP" => Ok(Self::ZP),
            _ => Ok(Self::Unknown),
        }
    }
}

impl std::fmt::Display for ControlledVocabulary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::ADO => "ADO",
                Self::AEO => "AEO",
                Self::AFO => "AFO",
                Self::AfPO => "AfPO",
                Self::AGRO => "AGRO",
                Self::AISM => "AISM",
                Self::AMPHX => "AMPHX",
                Self::APO => "APO",
                Self::APOLLO_SV => "APOLLO_SV",
                Self::ARO => "ARO",
                Self::BAO => "BAO",
                Self::BCGO => "BCGO",
                Self::BCO => "BCO",
                Self::BFO => "BFO",
                Self::BFO11 => "BFO11",
                Self::BIDO => "BIDO",
                Self::BMONT => "BMONT",
                Self::BSPO => "BSPO",
                Self::BTO => "BTO",
                Self::CARO => "CARO",
                Self::CCFO => "CCFO",
                Self::CDAO => "CDAO",
                Self::CDNO => "CDNO",
                Self::CEPH => "CEPH",
                Self::CHEBI => "CHEBI",
                Self::CHEMINF => "CHEMINF",
                Self::CHIRO => "CHIRO",
                Self::CHMO => "CHMO",
                Self::CIDO => "CIDO",
                Self::CIO => "CIO",
                Self::CL => "CL",
                Self::CLAO => "CLAO",
                Self::CLO => "CLO",
                Self::CLO_hPSCreg => "CLO-hPSCreg",
                Self::CLO_NICR => "CLO-NICR",
                Self::CLYH => "CLYH",
                Self::CMF => "CMF",
                Self::CMO => "CMO",
                Self::COB => "COB",
                Self::COLAO => "COLAO",
                Self::CRO => "CRO",
                Self::CTCAE_OAEview => "CTCAE-OAEview",
                Self::CTENO => "CTENO",
                Self::CTO => "CTO",
                Self::CVDO => "CVDO",
                Self::d_acts => "d-acts",
                Self::DDANAT => "DDANAT",
                Self::DDPHENO => "DDPHENO",
                Self::DIDEO => "DIDEO",
                Self::DINTO => "DINTO",
                Self::DISDRIV => "DISDRIV",
                Self::DOID => "DOID",
                Self::DRON => "DRON",
                Self::DUO => "DUO",
                Self::ECAO => "ECAO",
                Self::ECO => "ECO",
                Self::ECOCORE => "ECOCORE",
                Self::ECTO => "ECTO",
                Self::EDAM => "EDAM",
                Self::EFO => "EFO",
                Self::EHDAA2 => "EHDAA2",
                Self::EMAPA => "EMAPA",
                Self::ENVO => "ENVO",
                Self::EPIO => "EPIO",
                Self::EPO => "EPO",
                Self::ERO => "ERO",
                Self::EUPATH => "EUPATH",
                Self::EXMO => "EXMO",
                Self::ExO => "ExO",
                Self::FAO => "FAO",
                Self::FBbi => "FBbi",
                Self::FBbt => "FBbt",
                Self::FBcv => "FBcv",
                Self::FBdv => "FBdv",
                Self::FIDEO => "FIDEO",
                Self::FIX => "FIX",
                Self::FLOPO => "FLOPO",
                Self::FMA => "FMA",
                Self::FOBI => "FOBI",
                Self::FOODON => "FOODON",
                Self::FOVT => "FOVT",
                Self::FYPO => "FYPO",
                Self::GALLONT => "GALLONT",
                Self::GAZ => "GAZ",
                Self::GECKO => "GECKO",
                Self::GENEPIO => "GENEPIO",
                Self::GENO => "GENO",
                Self::GEO => "GEO",
                Self::GNO => "GNO",
                Self::GO => "GO",
                Self::GSSO => "GSSO",
                Self::HANCESTRO => "HANCESTRO",
                Self::HAO => "HAO",
                Self::HINO => "HINO",
                Self::HOM => "HOM",
                Self::HP => "HP",
                Self::HsapDv => "HsapDv",
                Self::HSO => "HSO",
                Self::HTN => "HTN",
                Self::IAO => "IAO",
                Self::IAO_Onto_Meta => "IAO-Onto-Meta",
                Self::ICDO => "ICDO",
                Self::ICEO => "ICEO",
                Self::ICO => "ICO",
                Self::IDO => "IDO",
                Self::IDOBRU => "IDOBRU",
                Self::IDOMAL => "IDOMAL",
                Self::INO => "INO",
                Self::KISAO => "KISAO",
                Self::KTAO => "KTAO",
                Self::LABO => "LABO",
                Self::LEPAO => "LEPAO",
                Self::LINCS_CLOview => "LINCS-CLOview",
                Self::LTHIDO => "LTHIDO",
                Self::MA => "MA",
                Self::MAMO => "MAMO",
                Self::MAXO => "MAXO",
                Self::MCO => "MCO",
                Self::MCRO => "MCRO",
                Self::MF => "MF",
                Self::MFMO => "MFMO",
                Self::MFOEM => "MFOEM",
                Self::MFOMD => "MFOMD",
                Self::MI => "MI",
                Self::MIAPA => "MIAPA",
                Self::MICRO => "MICRO",
                Self::MIDO => "MIDO",
                Self::miRNAO => "miRNAO",
                Self::MIRO => "MIRO",
                Self::MMO => "MMO",
                Self::MmusDv => "MmusDv",
                Self::MOD => "MOD",
                Self::MONDO => "MONDO",
                Self::MOP => "MOP",
                Self::MP => "MP",
                Self::MPATH => "MPATH",
                Self::MPIO => "MPIO",
                Self::MRO => "MRO",
                Self::MS => "MS",
                Self::NBO => "NBO",
                Self::NCBITaxon => "NCBITaxon",
                Self::NCIT => "NCIT",
                Self::NCRO => "NCRO",
                Self::NDF_RT => "NDF-RT",
                Self::NFO => "NFO",
                Self::NGBO => "NGBO",
                Self::NOMEN => "NOMEN",
                Self::OAE => "OAE",
                Self::OARCS => "OARCS",
                Self::OBA => "OBA",
                Self::OBCS => "OBCS",
                Self::OBI => "OBI",
                Self::OBI_NIAID_GSC_BRC_view => "OBI-NIAID-GSC-BRC-view",
                Self::OBIB => "OBIB",
                Self::OBIws => "OBIws",
                Self::OCCO => "OCCO",
                Self::OCE => "OCE",
                Self::OCMR => "OCMR",
                Self::OCRe => "OCRe",
                Self::OCVDAE => "OCVDAE",
                Self::ODAE => "ODAE",
                Self::ODNAE => "ODNAE",
                Self::OGG => "OGG",
                Self::OGG_At => "OGG-At",
                Self::OGG_Bru => "OGG-Bru",
                Self::OGG_Ce => "OGG-Ce",
                Self::OGG_Dm => "OGG-Dm",
                Self::OGG_Dr => "OGG-Dr",
                Self::OGG_Mm => "OGG-Mm",
                Self::OGG_Pf => "OGG-Pf",
                Self::OGG_Sc => "OGG-Sc",
                Self::OGI => "OGI",
                Self::OGMS => "OGMS",
                Self::OGSF => "OGSF",
                Self::OHD => "OHD",
                Self::OHMI => "OHMI",
                Self::OHPI => "OHPI",
                Self::OlatDv => "OlatDv",
                Self::OLOBO => "OLOBO",
                Self::OM => "OM",
                Self::OMIABIS => "OMIABIS",
                Self::OMIT => "OMIT",
                Self::OMO => "OMO",
                Self::OMP => "OMP",
                Self::OMRSE => "OMRSE",
                Self::ONE => "ONE",
                Self::ONS => "ONS",
                Self::ONTOAVIDA => "ONTOAVIDA",
                Self::ONTONEO => "ONTONEO",
                Self::OOSTT => "OOSTT",
                Self::OPL => "OPL",
                Self::OPMI => "OPMI",
                Self::ORNASEQ => "ORNASEQ",
                Self::OSCI => "OSCI",
                Self::OVAE => "OVAE",
                Self::PATO => "PATO",
                Self::PCL => "PCL",
                Self::PCO => "PCO",
                Self::PDRO => "PDRO",
                Self::PdumDv => "PdumDv",
                Self::PECO => "PECO",
                Self::PHASES => "PHASES",
                Self::PHIPO => "PHIPO",
                Self::PIDO => "PIDO",
                Self::PLANA => "PLANA",
                Self::PLANP => "PLANP",
                Self::PNO => "PNO",
                Self::PO => "PO",
                Self::PORO => "PORO",
                Self::PPO => "PPO",
                Self::PR => "PR",
                Self::PRIDE => "PRIDE",
                Self::PROCO => "PROCO",
                Self::PSDO => "PSDO",
                Self::PSO => "PSO",
                Self::PW => "PW",
                Self::RBO => "RBO",
                Self::REO => "REO",
                Self::RESID => "RESID",
                Self::REX => "REX",
                Self::RNAO => "RNAO",
                Self::RO => "RO",
                Self::RS => "RS",
                Self::RXNO => "RXNO",
                Self::SBO => "SBO",
                Self::SCDO => "SCDO",
                Self::SCO => "SCO",
                Self::SDGIO => "SDGIO",
                Self::SEP => "SEP",
                Self::SEPIO => "SEPIO",
                Self::SIBO => "SIBO",
                Self::SIO => "SIO",
                Self::SLSO => "SLSO",
                Self::SO => "SO",
                Self::SPD => "SPD",
                Self::STATO => "STATO",
                Self::SWO => "SWO",
                Self::SYMP => "SYMP",
                Self::T4FS => "T4FS",
                Self::TADS => "TADS",
                Self::TAXRANK => "TAXRANK",
                Self::TGMA => "TGMA",
                Self::TO => "TO",
                Self::TRANS => "TRANS",
                Self::TSO => "TSO",
                Self::TTO => "TTO",
                Self::TXPO => "TXPO",
                Self::UBERON => "UBERON",
                Self::uberon_basic => "uberon-basic",
                Self::uberon_collected_metazoa => "uberon_collected_metazoa",
                Self::uberon_depictions => "uberon_depictions",
                Self::UO => "UO",
                Self::UPA => "UPA",
                Self::UPHENO => "UPHENO",
                Self::VariO => "VariO",
                Self::VBO => "VBO",
                Self::VICO => "VICO",
                Self::VIDO => "VIDO",
                Self::VIO => "VIO",
                Self::VIVO_ISF => "VIVO-ISF",
                Self::VO => "VO",
                Self::VT => "VT",
                Self::VTO => "VTO",
                Self::WBbt => "WBbt",
                Self::WBls => "WBls",
                Self::WBPhenotype => "WBPhenotype",
                Self::XAO => "XAO",
                Self::XCO => "XCO",
                Self::XL => "XL",
                Self::XLMOD => "XLMOD",
                Self::XPO => "XPO",
                Self::ZECO => "ZECO",
                Self::ZFA => "ZFA",
                Self::ZFS => "ZFS",
                Self::ZP => "ZP",
                Self::Unknown => "?",
            }
        )
    }
}

/// A CURIE is a namespace + accession identifier
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Curie {
    /// The controlled vocabulary
    pub cv: ControlledVocabulary,
    /// The accession code
    pub accession: AccessionCode,
}

/// A curie. For documentation reasons a full term definition is also allowed. So both examples
/// give the exact same `Curie` result:
/// ```
/// use mzcv::curie;
/// curie!(MS:1002357);
/// curie!(MS:1002357|PSM-level probability);
/// ```
#[macro_export]
macro_rules! curie {
    ($ns:ident:$acc:tt) => {
        $crate::Curie {
            cv: $crate::ControlledVocabulary::$ns,
            accession: $crate::accession_code!($acc),
        }
    };
    ($ns:ident:$acc:literal|$($name:tt)+) => {
        $crate::curie!($ns:$acc)
    };
}

impl std::fmt::Display for Curie {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}", self.cv, self.accession)
    }
}

/// An error that occured when parsing a CURIE
#[derive(Debug)]
pub enum CURIEParsingError {
    /// The CV name was not recognised
    UnknownControlledVocabulary,
    /// The accession code could not be parsed
    AccessionParsingError(AccessionCodeParseError),
    /// No name separator ':' could be found
    MissingNamespaceSeparator,
}

impl std::str::FromStr for Curie {
    type Err = CURIEParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((cv, accession)) = s.split_once(':').or_else(|| s.split_once('_')) {
            let cv = cv.parse::<ControlledVocabulary>().unwrap(); // Unknown CVs are handled with the CV::Unknown option
            let accession = accession
                .parse()
                .map_err(CURIEParsingError::AccessionParsingError)?;
            Ok(Self { cv, accession })
        } else {
            Err(CURIEParsingError::MissingNamespaceSeparator)
        }
    }
}

/// An accession code, Can either be a numeric code (u32 to 4 milion, so 9 fully utilised digits).
/// Or it can be an ASCII alphanumeric code (case-sensitive) of 1 to 8 characters.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AccessionCode {
    /// A strictly numeric code
    Numeric(u32),
    /// An alphanumeric code, can be length 1 to 8 ASCII chars (no spaces allowed)
    Alphanumeric(NonZeroU8, [u8; 7]),
    // Maybe use ThinStr for dynamic strings, very likely this cannot be done together with an 8 byte inline string but it might be better in most ways
    //Big(ThinStr),
}

impl std::fmt::Display for AccessionCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Numeric(n) => write!(f, "{n:07}"),
            Self::Alphanumeric(first, bytes) => {
                write!(
                    f,
                    "{}{}",
                    char::from(first.get()),
                    std::str::from_utf8(bytes.trim_ascii_end()).unwrap()
                )
            }
        }
    }
}

/// Create an accession code, using this with a numeric code is match arm valid, using this with an alphanumeric code is const compatible.
/// ```rust ignore
/// match curie!(MS:000111) {
///     curie!(MS:0000111) => (),
///     x if x == curie!(MS:H000111G) => (),
///     _ => (),
/// }
/// ```
#[macro_export]
macro_rules! accession_code {
    ($acc:literal) => {
        // #[allow(clippy::zero_prefixed_literal)] cannot put allows on single expressions yet
        $crate::AccessionCode::Numeric($acc)
    };
    ($acc:tt) => {
        const {
            const BYTES: &[u8] = stringify!($acc).as_bytes();
            assert!(BYTES[0] != 0);
            assert!(BYTES.is_ascii());
            match BYTES {
                [a, b, c, d, e, f, g, h] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, *d, *e, *f, *g, *h],
                ),
                [a, b, c, d, e, f, g] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, *d, *e, *f, *g, 0x20],
                ),
                [a, b, c, d, e, f] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, *d, *e, *f, 0x20, 0x20],
                ),
                [a, b, c, d, e] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, *d, *e, 0x20, 0x20, 0x20],
                ),
                [a, b, c, d] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, *d, 0x20, 0x20, 0x20, 0x20],
                ),
                [a, b, c] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, 0x20, 0x20, 0x20, 0x20, 0x20],
                ),
                [a, b] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20],
                ),
                [a] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20],
                ),
                _ => panic!(concat!(
                    "Cannot convert ",
                    stringify!($acc),
                    " to accession code. Can only handle strings of 1 to 8 bytes"
                )),
            }
        }
    };
}

/// An error when parsing an accession code
#[derive(Debug)]
pub enum AccessionCodeParseError {
    /// An invalid character was encountered, only alphanumeric characters are allowed
    InvalidCharacters(String),
    /// The accession code was over 8 bytes
    TooLong(String),
    /// The accession code was empty
    Empty,
}

impl std::str::FromStr for AccessionCode {
    type Err = AccessionCodeParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(AccessionCodeParseError::Empty);
        }
        if s.chars().any(|c| !c.is_alphanumeric()) {
            return Err(AccessionCodeParseError::InvalidCharacters(s.to_string()));
        }
        s.parse::<u32>().map_or_else(
            |_| {
                if s.len() > 8 {
                    return Err(AccessionCodeParseError::TooLong(s.to_string()));
                }
                let mut bytes = [0x20; 7]; // Space
                bytes[..s.len() - 1].copy_from_slice(&s.as_bytes()[1..]);
                curie!(NCIT:P383);
                Ok(Self::Alphanumeric(
                    NonZeroU8::new(s.as_bytes()[0]).unwrap(),
                    bytes,
                ))
            },
            |u| Ok(Self::Numeric(u)),
        )
    }
}

/// A term, a CURIE plus its name
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Term {
    /// The CURIE (eg MS:0000000)
    pub accession: Curie,
    /// The name
    pub name: std::borrow::Cow<'static, str>, // Static Cow to allow to have term in match arms
}

/// Create a new term `term!(MS:1002357|PSM-level probability)`. The accession/name combination
/// is not validated.
#[macro_export]
macro_rules! term {
    ($ns:ident:$accession:literal|$($name:tt)+) =>  {
        $crate::mzspeclib::Term {
            accession: $crate::curie!($ns:$accession),
            name: std::borrow::Cow::Borrowed(stringify!($($name)+))
        }
    };
}

#[cfg(test)]
#[allow(
    clippy::missing_panics_doc,
    clippy::zero_prefixed_literal,
    clippy::unnested_or_patterns
)]
mod tests {
    use crate::Curie;

    #[test]
    fn parse_curies() {
        let options = [
            "NCIT:C25330",
            "NCIT:R100",
            "NCIT:P383",
            "MS:1000014",
            "UO:0000245",
            "BFO:0000015",
            "BAO_0000925", // Seen in one Obo file, but also quite common in URLs
            "GAZ:00002933",
            "AfPO_0000233",
            "BTO:0004947",
            "PRIDE:0000521",
            "MOD:01188",
            "XLMOD:07097",
            "GNO:G00001NT",
            "GNO:00000015",
        ];
        for option in options {
            let _curie: Curie = option.parse().unwrap();
        }
    }

    #[test]
    fn curie_macro() {
        assert_eq!("NCIT:C25330".parse::<Curie>().unwrap(), curie!(NCIT:C25330));
        assert_eq!("NCIT:R100".parse::<Curie>().unwrap(), curie!(NCIT:R100));
        assert_eq!("NCIT:P383".parse::<Curie>().unwrap(), curie!(NCIT:P383));
        assert_eq!("MS:1000014".parse::<Curie>().unwrap(), curie!(MS:1000014));
        assert_eq!("UO:0000245".parse::<Curie>().unwrap(), curie!(UO:0000245));
        assert_eq!("BFO:0000015".parse::<Curie>().unwrap(), curie!(BFO:0000015));
        assert_eq!("BTO:0004947".parse::<Curie>().unwrap(), curie!(BTO:0004947));
        assert_eq!(
            "PRIDE:0000521".parse::<Curie>().unwrap(),
            curie!(PRIDE:0000521)
        );
        assert_eq!("MOD:01188".parse::<Curie>().unwrap(), curie!(MOD:01188));
        assert_eq!("XLMOD:07097".parse::<Curie>().unwrap(), curie!(XLMOD:07097));
        assert_eq!(
            "GNO:G00001NT".parse::<Curie>().unwrap(),
            curie!(GNO:G00001NT)
        );
        assert_eq!(
            "GNO:00000015".parse::<Curie>().unwrap(),
            curie!(GNO:00000015)
        );
    }

    #[test]
    fn match_arm() {
        assert!(
            match "NCIT:C25330".parse::<Curie>().unwrap() {
                curie!(MS:0000111) | curie!(MS:0000112|hello world) => Err(()),
                x if x == curie!(MS:H000111G) => Err(()),
                x if x == curie!(NCIT:C25330) => Ok(()),
                _ => Err(()),
            }
            .is_ok()
        );
    }
}
