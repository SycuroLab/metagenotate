# ***************************************
# * Snakefile for met assemble pipeline *
# ***************************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
import os

# Set the PATH environment variable for metaWRAP bin folder.
os.environ["PATH"]+=os.pathsep+"/bulk/IMCshared_bulk/shared/shared_software/metaWRAP/bin"

SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input: 
        expand(config["output_dir"]+"/{sample}/assembly/final_assembly.fasta",sample=SAMPLES)

rule metawrap_assembly:
    input:
        fastq_read1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        fastq_read2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly"),
        metagenome_final_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","final_assembly.fasta")
    params: 
        memory_in_gb = config["memory_in_gb"],
        threads = config["assembler_threads"],
        min_contig_length = config["min_contig_length"]
    conda: "utils/envs/metawrap_env.yaml"
    shell:
        "metawrap assembly -1 {input.fastq_read1} -2 {input.fastq_read2} -m {params.memory_in_gb} -t {params.threads} -l {params.min_contig_length} --metaspades -o {output.sample_assembly_dir}"


rule metawrap_binning:
    input:
         metagenome_final_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","final_assembly.fasta"),
         fastq_read12 = os.path.join(config["input_dir"],"{sample}*fastq")
    output:
         metabat2_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2_bins"),
         maxbin2_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2_bins"),
         concoct_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct_bins")

    params:
         sample_initial_binning_dir = os.path.join(config["output_dir"],"{sample}","initial_binning"),
         threads = config["binning_threads"],
         renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
    conda: "utils/envs/metawrap_env.yaml"
    shell:
         "cp {input.metagenome_final_assembly_file} {params.renamed_metagenome_assembly_file} "
         "metawrap binning -o {params.sample_initial_binning_dir} -t {params.threads} -a {params.renamed_metagenome_assembly_file} --metabat2 --maxbin2 --concoct {input.fastq_read12}"

