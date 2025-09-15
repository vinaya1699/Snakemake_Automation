# Snakemake Automation
This is an automated snakemake pipeline for transcriptomics analysis . Currently , this script works till read count generation.

# To install snakemake 
conda create -c conda-forge -c bioconda -n snakemake snakemake

# To check the succesful jobs
snakemake -n

# To run snakemake pipeline
snakemake -j 2

Options:
-j : Parallel Jobs
