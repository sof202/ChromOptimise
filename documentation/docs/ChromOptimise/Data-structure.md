---
sidebar_position: 2
---

# Data structure

A guide for the structure of the data directory is given below (You only need 
to create the directories starting with an integer):

```text
Main_Data_Directory
├── YourProcessedBamAndBedFiles
│   ├── Epigenetic_Mark_1.bed
│   ├── Epigenetic_Mark_2.bam
│   ├── ...
│   └── Epigenetic_Mark_n.bed
├── 1_Binary_Files
│   └── Results_From_Run_n
│       └── Binary files for each chromosome
├── 2_Model_Files
│   └── Results_From_Run_n
│       ├── Likelihood_Values
│       │   └── likelihoods.txt 
│       ├── STATEBYLINE
│       │   └── cell_type.statebyline.txt 
│       └── ChromHMM model files
├── 3_Optimum_Number_Of_States
│   └── Results_From_Run_n
│       ├── Euclidean_distances
│       │   └── Euclidean_distances_model-n.txt 
│       ├── Flanking_states
│       │   └── Likeliest_flanking_states_model-n.txt 
│       ├── Isolation_scores
│       │   └── Isolation_scores_model-n.txt 
│       ├── Redundant_states_model-n.txt 
│       └── OptimumNumberOfStates.txt
├── 4_LDSC_Assessment_Files
│   └── Results_From_Run_n
│       ├── annotation
│       │   ├── ChromHMM.n.l2.ldscore.gz
│       │   ├── ChromHMM.n.annot
│       │   ├── ChromHMM.n.log
│       │   ├── ChromHMM.n.M
│       │   └── ChromHMM.n.M_5_50
│       ├── heritability
│       │   ├── gwas_trait_n.results
│       │   └── gwas_trait_n.log
│       └── plots
│           ├── All_categories
│           │   ├── Enrichment_heatmap.png
│           │   └── Enrichment_pvalues_gwas_trait_n.png
│           └── State_categories
│               ├── Enrichment_heatmap.png
│               └── Enrichment_pvalues_gwas_trait_n.png
└── 5_Big_Model_Files
    ├── ChromHMM model files
    └── Plots
        ├── Euclidean_Distance_Histrograms
        └── Transition_Maxima_Scatter_Plots


LDSC_reference_files
├── PLINK_files
├── frq_files 
├── weights_files
├── baseline_files 
└── gwas_traits
```

