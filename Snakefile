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
        expand(config["output_dir"]+"/{sample}/assembly/metaspades/scaffolds.fasta",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/final_assembly.fasta",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/{sample}_metagenome.fasta",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/initial_binning/metabat2_bins/bin.1.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/initial_binning/maxbin2_bins/bin.0.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/initial_binning/concoct_bins/bin.0.fa",sample=SAMPLES),

rule metawrap_assembly:
    input:
        fastq_read1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        fastq_read2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        metagenome_final_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","final_assembly.fasta"),
        metaspades_scaffolds_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/scaffolds.fasta"),
        assembly_report_file = os.path.join(config["output_dir"],"{sample}","assembly","assembly_report.html")
    params: 
        memory_in_gb = config["memory_in_gb"],
        threads = config["assembler_threads"],
        min_contig_length = config["min_contig_length"],
        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly")
    conda: "utils/envs/metawrap_env.yaml"
    shell:
        "metawrap assembly -1 {input.fastq_read1} -2 {input.fastq_read2} -m {params.memory_in_gb} -t {params.threads} -l {params.min_contig_length} --metaspades -o {params.sample_assembly_dir}"

rule rename_metawrap_assembly_file:
    input:
        metagenome_final_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","final_assembly.fasta")    
    output:
        renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
    params:
        threads = config["binning_threads"],
        fastq_read12 = os.path.join(config["input_dir"],"{sample}*fastq")
    conda: "utils/envs/metawrap_env.yaml"
    shell:
        "cp {input.metagenome_final_assembly_file} {output.renamed_metagenome_assembly_file}"


rule metawrap_binning:
    input:
         renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
    output:
         metabat2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2_bins","bin.1.fa"),
         maxbin2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2_bins","bin.0.fa"),
         concoct_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct_bins","bin.0.fa")
    params:
         sample_initial_binning_dir = os.path.join(config["output_dir"],"{sample}","initial_binning"),
         threads = config["binning_threads"],
         fastq_read12 = os.path.join(config["input_dir"],"{sample}*fastq")
    conda: "utils/envs/metawrap_env.yaml"
    shell:
         "metawrap binning -o {params.sample_initial_binning_dir} -t {params.threads} -a {input.renamed_metagenome_assembly_file} --metabat2 --maxbin2 --concoct {params.fastq_read12}"

rule metawrap_bin_refinement:
    input:
         metabat2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2_bins","bin.1.fa"),
         maxbin2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2_bins","bin.0.fa"),
         concoct_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct_bins","bin.0.fa")

    output:
         sample_bin_refinement_dir = os.path.join(config["output_dir"],"{sample}","initial_binning"),
    params:
         threads = config["bin_refinement_threads"],
         metabat2_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2_bins"),
         maxbin2_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2_bins"),
         concoct_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct_bins"),
         completeness_thresh = config["completeness_thresh"],
         contamination_thresh = config["contamination_thresh"]
    conda: "utils/envs/metawrap_env.yaml"
    shell:
         "metawrap bin_refinement -o {output.sample_bin_refinement_dir} -t {params.threads} -A {params.metabat2_bins_dir} -B {params.maxbin2_bins_dir} -C {params.concoct_bins_dir} -c {params.completeness_thresh} -x {params.contamination_thresh}"

