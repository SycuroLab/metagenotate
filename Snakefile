# ***************************************
# * Snakefile for met assemble pipeline *
# ***************************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
import os

SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input:
        os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
        

rule merge_reads:
    input:
        fastq_read1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        fastq_read2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly")
    params: 
        metaphlan_database = config["metaphlan_database"],
        memory_in_gb = config["memory_in_gb"]
        threads = config["assembler_threads"],
        min_contig_length = config["min_contig_length"]
    conda: "utils/envs/metawrap_env.yaml"
    shell:
        "metawrap assembly -1 {input.fastq_read1} -2 {input.fastq_read2} -m {params.memory_in_gb} -t {params.threads} -l {params.min_contig_length} --metaspades -o {output.sample_assembly_dir}"
