# Snakemake Automation
This is an automated snakemake pipeline for transcriptomics analysis . Currently , this script works till read count generation.

# To install snakemake 
conda create -c conda-forge -c bioconda -n snakemake snakemake

# To use this snakemake script
git clone https://github.com/vinaya1699/Snakemake_Automation.git

# To check the succesful jobs
snakemake -n

# To run snakemake pipeline
snakemake -j 2

Options:
-j : Parallel Jobs

# Example Config File :
samples:
  Zinc_Sample: Zinc_Sample
  Control_Sample: Control_Sample

organism:
  Mus_musculus: Mus_musculus
  
threads:
  threads: 10
