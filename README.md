## Introduction
This is a standard GWAS pre-processing pipeline which includes GWAS, meta-analysis, genetic correlation/heritability, Mendelian Randomisation, polygenic risk score prediction and gene prioritisation.  
Everything is processed on the basis of **groups of phenotypes**.

## Usage
1. The starting point is `pheno_batch.py` which batch runs a `pheno.py` to be defined separately.
2. This `pheno_batch.py` script generates a per-subject phenotype file in tabular format and can be concatenated using `pheno_concat.py`. This results in a TXT file in the following format: `FID  IID  <pheno_0>  <pheno_1> ...`. All phenotypes in this file are considered a group of phenotype. The name of this file should be used for the rest of the pipeline.
3. Call `gwa_batch.py <name of phenotype group>` to run fastGWA-MLM for every phenotype in that group.
4. External GWAS summary statistics can be integrated into the pipeline at this stage; they should be re-formatted to the fastGWA format with columns: `CHR SNP POS A1 A2 AF1 N BETA/OR SE P` and optional `INFO N_CAS N_CON DIRE` columns. A `metadata` file can be manually included in every phenotype group folder for information that cannot be directly inferred from the GWAS summary statistics (e.g. effective sample size for genomic SEM).
5. For all other scripts, call `*.py --help` to get complete usage options and parameters. For all `pheno` arguments, use **groups of phenotypes** separated by blank spaces.
6. All scripts with the `_batch` label require a SLURM manager so that lengthy jobs can be run parallel in the background.
7. All scripts with the `_parse` label generate human-readable summary tables.

## Order of scripts
1. `pheno_*.py`
2. `gwa_*.py`
3. Correlative analysis pipeline
    1. `heri_batch.py`: Heritability analysis using LDSC
    2. `gcorr_batch.py`: Genetic correlation analysis using LDSC
    3. `gwa_clump_batch.py`: Clump for independent loci at user-defined significance thresholds
    4. `mr_*.py`: Mendelian Randomisation analysis
        1. `mr_extract_snp_batch.py`
        2. `mr_batch.py`
        3. `mr_mvmr_batch.py`
4. Gene annotation pipeline
    1. `annot_magma_batch.py`: Annotate using MAGMA and H-MAGMA to aggregate SNP to gene-level summary stats
    2. `annot_smr_batch.py`: Annotate using summary data Mendelian Randomisation to combine GWAS summary stats with xQTL data
    3. `finemap_*.py`: Fine-map for causal variants and then annotate using Ensembl-VEP
5. Polygenic score pipeline
    1. `prs_from_gwa_batch.py`, `prs_score_batch.py`: Generate SNP-wise effect sizes using PRScs and then score for the target cohort
    2. `prs_corr_*.py`: Correlation analysis
6. Genomic SEM pipeline
    1. `gsem_batch.py` is a very flexible framework which can implement common factor analysis, EFA/CFA, causal models, GWAS by subtraction and manually input models.