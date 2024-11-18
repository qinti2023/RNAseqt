import os
CONFIG_YAML = """\
samples:
  sample1: 
    R1: "path to sample1 R1" 
    R2: "path to sample1 R2" 
  sample2: 
    R1: "path to sample2 R1" 
    R2: "path to sample2 R2" 
reference: "Reference/mm10/genome"  
gtf: "Reference/gencode.vM25.primary_assembly.annotation.gtf.gz"  
counts: "Counts"  
threads: 8   
"""
SNAKEFILE = """\
import glob
configfile: "config.yaml"
rule all:
    input:
        "QCresult/multiqc_report.html",  
        expand("Cleandata/{sample}_R1_val_1.fq", sample=config["samples"]),
        expand("Cleandata/{sample}_R2_val_2.fq", sample=config["samples"]),
        "Mapping/multiqc_report.html",  
        "Counts/counts.txt",             
        "Cleandata/multiqc_report.html"  
rule rename:
    input:
        r1=lambda wildcards: config["samples"][wildcards.sample]["R1"],
        r2=lambda wildcards: config["samples"][wildcards.sample]["R2"]
    output:
        r1_new = "Rawdata/{sample}_R1.fq",
        r2_new = "Rawdata/{sample}_R2.fq"
    shell:
        "mv {input.r1} {output.r1_new} && mv {input.r2} {output.r2_new}"
rule fastqc:
    input:
        r1_new = "Rawdata/{sample}_R1.fq",
        r2_new = "Rawdata/{sample}_R2.fq"
    output:
        r1_zip = "QCresult/{sample}_R1_fastqc.zip",
        r1_html = "QCresult/{sample}_R1_fastqc.html",
        r2_zip = "QCresult/{sample}_R2_fastqc.zip",
        r2_html = "QCresult/{sample}_R2_fastqc.html"
    threads: config["threads"]
    shell:
        "fastqc -t {threads} -o QCresult {input.r1_new} {input.r2_new}"

rule qc_multiqc:
    input:
        html=expand("QCresult/{sample}_R1_fastqc.html", sample=config["samples"]) +
             expand("QCresult/{sample}_R2_fastqc.html", sample=config["samples"]),
        zip=expand("QCresult/{sample}_R1_fastqc.zip", sample=config["samples"]) +
            expand("QCresult/{sample}_R2_fastqc.zip", sample=config["samples"])
    output:
        report="QCresult/multiqc_report.html"
    shell:
        "multiqc -o QCresult QCresult/*_fastqc.zip"
        
rule trim_galore:
    input:
        r1="Rawdata/{sample}_R1.fq",
        r2="Rawdata/{sample}_R2.fq"
    output:
        r1_clean="Cleandata/{sample}_R1_val_1.fq",
        r2_clean="Cleandata/{sample}_R2_val_2.fq",
        r1_clean_zip = "Cleandata/{sample}_R1_val_1_fastqc.zip",
        r1_clean_html = "Cleandata/{sample}_R1_val_1_fastqc.html",
        r2_clean_zip = "Cleandata/{sample}_R2_val_2_fastqc.zip",
        r2_clean_html = "Cleandata/{sample}_R2_val_2_fastqc.html"
    threads: config["threads"]
    shell:
        "trim_galore -q 20 --length 36 --max_n 3 --stringency 3 --fastqc --paired -o Cleandata {input.r1} {input.r2}"

rule clean_multiqc:
    input:
        html=expand("Cleandata/{sample}_R1_val_1_fastqc.html", sample=config["samples"]) +
             expand("Cleandata/{sample}_R2_val_2_fastqc.html", sample=config["samples"]),
        zip=expand("Cleandata/{sample}_R1_val_1_fastqc.zip", sample=config["samples"]) +
            expand("Cleandata/{sample}_R2_val_2_fastqc.html", sample=config["samples"])
    output:
        report="Cleandata/multiqc_report.html"
    shell:
        "multiqc -o Cleandata Cleandata/*_fastqc.zip"
        
rule hisat2:
    input:
        r1_clean="Cleandata/{sample}_R1_val_1.fq",
        r2_clean="Cleandata/{sample}_R2_val_2.fq"
    output:
        bam="Mapping/{sample}.Hisat_aln.sorted.bam",
        bai="Mapping/{sample}.Hisat_aln.sorted.bam.bai"
    log:
        "Mapping/{sample}.log"
    threads: config["threads"]
    shell:
        "hisat2 -p {threads} -x {config[reference]} -1 {input.r1_clean} -2 {input.r2_clean} 2> {log} | samtools sort -@ {threads} -o {output.bam} &&"
        "samtools index {output.bam}"
        
rule map_multilog:
    input:
        log=expand("Mapping/{sample}.log", sample=config["samples"])
    output:
        report="Mapping/multiqc_report.html"
    shell:
        "multiqc -o Mapping {input.log}"  

rule featureCounts:
    input:
        bams=expand("Mapping/{sample}.Hisat_aln.sorted.bam", sample=config["samples"]),
        gtf=config["gtf"]
    output:
        counts="Counts/all.id.txt",
        summary="Counts/all.id.txt.summary"
    threads: config["threads"]
    shell:
        "featureCounts -T {threads} -p -t exon -g gene_id -a {input.gtf} -o {output.counts} {input.bams}"
rule counts_matrix:
    input:
        "Counts/all.id.txt"
    output:
        "Counts/counts.txt"
    shell:
        "cat {input} | cut -f1,7- > {output}"
"""
def init_project(base_dir):
    """
    Create necessary files and folders for the RNA-seq pipeline.
    """
    base_dir = os.path.abspath(base_dir)
    folders = ["Rawdata", "Cleandata", "QCresult", "Mapping", "Counts", "Reference"]
    for folder in folders:
        os.makedirs(os.path.join(base_dir, folder), exist_ok=True)
    with open(os.path.join(base_dir, "config.yaml"), "w") as f:
        f.write(CONFIG_YAML)
    with open(os.path.join(base_dir, "Snakefile"), "w") as f:
        f.write(SNAKEFILE)
    print(f"Pipeline initialized in: {base_dir}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate Snakemake project structure for RNA-seq.")
    parser.add_argument("base_dir", help="The root directory where the pipeline will be initialized.")
    args = parser.parse_args()
    init_project(args.base_dir)


