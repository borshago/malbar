Contents:
 - Snakefile4: the snakemake workflow script
 - config4.yml: the configuration file for Snakefile4
 - malbar1.yml: the specifications for creating the conda environment malbar1
 - malbar1_export.yml: the export of the full malbar1 conda environment
 - malbar2.yml: the specifications for creating the conda environment malbar2
 - malbar2_export.yml: the export of the full malbar2 conda environment
 - CODE/analysis.r: the R analysis of tag and mRNA counts
 - CODE/demux_seurat.r: the R analysis of tags with Seurat HTODemux()


The scripts expect the following directory structure:

```
CODE
DATA
├── CHROMIUM
│   ├── P37252_1011
│   ├── P37252_1021
│   ├── X3_24_038_GE
│   └── X3_24_039_GE
├── FASTQ
└── TRIMMED
INTERMEDIATE
├── CITESEQ_COUNT
│   ├── P32860_1002_S3_L001
│   │   ├── EMPTYDROPS
│   ├── P32860_1004_S5_L001
│   │   ├── EMPTYDROPS
│   ├── P37252_1015_S2_L005
│   │   ├── EMPTYDROPS
│   └── P37252_1025_S4_L005
│       ├── EMPTYDROPS
├── FASTQC
├── R_ANALYSIS
└── SEURAT
LOGS
RESULTS
Snakefile4
config4.yml
```
