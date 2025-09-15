from datetime import datetime

print(f"\n🚀 Here begins the automated transcriptomics analysis pipeline! [{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\n")

configfile: "config.yaml"

SAMPLES = list(config["samples"].keys())
ORGANISM = list(config["organism"].keys())

rule all:
    input:
        # Reference index
        expand("0_Reference_Genome/{organism}.{i}.ht2", organism=ORGANISM, i=range(1,9)),
        
        # Raw FastQC outputs
        expand("1_RawData_Fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("1_RawData_Fastqc/{sample}_R2_fastqc.html", sample=SAMPLES),

        # Cleaned FASTQ files
        expand("2_CleanData/{sample}_R1.fq.gz", sample=SAMPLES),
        expand("2_CleanData/{sample}_R2.fq.gz", sample=SAMPLES),

        # Cleaned FastQC outputs
        expand("2_CleanData_Fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("2_CleanData_Fastqc/{sample}_R2_fastqc.html", sample=SAMPLES),
        
        # Alignment outputs
        expand("3_Alignment/{sample}.bam", sample=SAMPLES),
        
        # RC Input outputs
        expand("3_Alignment/Bam_files.txt"),
        
        # Read Count Generation
        expand("4_Featurecount_WTA/counts.txt")
        
rule Indexing_of_ref_Genome:
    input:
        "0_Reference_Genome/{organism}.fasta"
    output:
        expand("0_Reference_Genome/{{organism}}.{i}.ht2", i=range(1,9))
    threads: 8
    shell:
        "hisat2-build {input} 0_Reference_Genome/{wildcards.organism}"

rule Quality_check:
    input:
        r1 = "1_RawData/{sample}_R1.fastq.gz",
        r2 = "1_RawData/{sample}_R2.fastq.gz"
    output:
        html_r1 = "1_RawData_Fastqc/{sample}_R1_fastqc.html",
        html_r2 = "1_RawData_Fastqc/{sample}_R2_fastqc.html",
        zip_r1  = "1_RawData_Fastqc/{sample}_R1_fastqc.zip",
        zip_r2  = "1_RawData_Fastqc/{sample}_R2_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    threads: 10
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} -o 1_RawData_Fastqc/ > {log} 2>&1
        """

rule Trimming:
    input:
        r1 = "1_RawData/{sample}_R1.fastq.gz",
        r2 = "1_RawData/{sample}_R2.fastq.gz"
    output:
        r1_clean = "2_CleanData/{sample}_R1.fq.gz",
        r2_clean = "2_CleanData/{sample}_R2.fq.gz",
        html = "2_CleanData/{sample}.html",
        json = "2_CleanData/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    threads: 20
    shell:
        """
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1_clean} \
            -O {output.r2_clean} \
            -h {output.html} \
            -j {output.json} \
            -w {threads} \
            -c \
            --detect_adapter_for_pe \
            > {log} 2>&1
        """

        
rule Quality_check_clean_files:
    input:
        r1 = "2_CleanData/{sample}_R1.fq.gz",
        r2 = "2_CleanData/{sample}_R2.fq.gz"
    output:
        html_r1 = "2_CleanData_Fastqc/{sample}_R1_fastqc.html",
        html_r2 = "2_CleanData_Fastqc/{sample}_R2_fastqc.html",
        zip_r1  = "2_CleanData_Fastqc/{sample}_R1_fastqc.zip",
        zip_r2  = "2_CleanData_Fastqc/{sample}_R2_fastqc.zip"
    log:
        "logs/fastqc/{sample}_clean.log"
    threads: 10
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} -o 2_CleanData_Fastqc/ > {log} 2>&1
        """
        
rule Alignment_of_clean_files:
    input:
        r1 = "2_CleanData/{sample}_R1.fq.gz",
        r2 = "2_CleanData/{sample}_R2.fq.gz"
    output:
        mapped = "3_Alignment/{sample}.bam",
        stats = "3_Alignment/{sample}_Mapping_Stat.txt"
        
    log:
        "logs/Alignment/{sample}_Alignment.log"
    threads: 10
    shell:
        """
        hisat2 -p {threads} -x 0_Reference_Genome/{ORGANISM} -1 {input.r1} -2 {input.r2} | samtools sort -@ {threads} -o {output.mapped} > {log} 2>&1
        samtools flagstat {output.mapped} > {output.stats}
        """

rule Raw_Read_Count_Input:
    input:
        bams = expand("3_Alignment/{sample}.bam", sample=SAMPLES)
    output:
        sample_names = "3_Alignment/Bam_files.txt"
    threads: 1
    shell:
        """
        ls {input.bams} > {output.sample_names}
        """
        
rule Raw_Read_Count_Generation:
    input:
        input_RC = "3_Alignment/Bam_files.txt"
    output:
        RC = "4_Featurecount_WTA/counts.txt"
    log:
        "logs/Featurecount/Featurecount.log"
    threads: 10
    params:
        gtf = lambda wildcards: f"0_Reference_Genome/{list(config['organism'].keys())[0]}.gtf"
    shell:
        """
        featureCounts -T {threads} -p -t gene -g gene_id -a {params.gtf} -o {output.RC} $(cat {input.input_RC}) > {log} 2>&1
        """
        
print(f"\n🚀 Here ends the automated transcriptomics analysis pipeline! [{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\n")
