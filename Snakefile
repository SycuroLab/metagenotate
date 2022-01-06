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
        expand(config["output_dir"]+"/{sample}/assembly/quast/transposed_report.tsv",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/prokka/prokka_annotations/{sample}_metagenome.fna",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/prokka/prokka_annotations/{sample}_metagenome.gff",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/metaerg/data/all.gff",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/initial_binning/metabat2_bins/bin.1.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/initial_binning/maxbin2_bins/bin.0.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/initial_binning/concoct_bins/bin.0.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/bin_refinement/metabat2_bins/bin.1.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/bin_refinement/maxbin2_bins/bin.0.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/bin_refinement/concoct_bins/bin.0.fa",sample=SAMPLES)

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

rule quast_assembly:
    input:
        renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
    output:
        quast_transposed_report_file = os.path.join(config["output_dir"],"{sample}","assembly","quast","transposed_report.tsv")
    params:
        assembly_quast_dir = os.path.join(config["output_dir"],"{sample}","assembly","quast"),
        threads = config["quast_threads"]
    conda: "utils/envs/quast_env.yaml"
    shell:
       "quast.py --output-dir {params.assembly_quast_dir} --threads {params.threads} {input.renamed_metagenome_assembly_file}"

rule prokka_assembly:
    input:
        renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
    output:
        prokka_fna_file = os.path.join(config["output_dir"],"{sample}","assembly","prokka","prokka_annotations","{sample}_metagenome.fna"),
        prokka_gff_file = os.path.join(config["output_dir"],"{sample}","assembly","prokka","prokka_annotations","{sample}_metagenome.gff")
    params:
        assembly_prokka_dir = os.path.join(config["output_dir"],"{sample}","assembly","prokka"),
        threads = config["prokka_threads"]
    conda: "utils/envs/prokka_env.yaml"
    shell:
       "python utils/scripts/prokka.py --genome_fasta_infile {input.renamed_metagenome_assembly_file} --num_cpus {params.threads} --metagenome true --output_dir {params.assembly_prokka_dir}"

rule metaerg_assembly:
    input:
        renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
    output:
#        metaerg_fna_file = os.path.join(config["output_dir"],"{sample}","assembly","metaerg","{sample}_metagenome.fna"),
        metaerg_gff_file = os.path.join(config["output_dir"],"{sample}","assembly","metaerg","data","all.gff")
    params:
        assembly_metaerg_dir = os.path.join(config["output_dir"],"{sample}","assembly","metaerg"),
        metaerg_database_path = config["metaerg_database_path"],
        locustag = "{sample}_metagenome",
        threads = config["metaerg_threads"]
    shell:
       "singularity run -H $HOME -B {params.metaerg_database_path}:/NGStools/metaerg/db -B /work:/work -B /bulk:/bulk /global/software/singularity/images/software/metaerg2.sif /NGStools/metaerg/bin/metaerg.pl --mincontiglen 200 --gcode 11 --gtype meta --minorflen 180 --cpus {params.threads} --evalue 1e-05 --identity 20 --coverage 70 --locustag {params.locustag} --force --outdir {params.assembly_metaerg_dir} {input.renamed_metagenome_assembly_file}"
 
#rule metawrap_binning:
#    input:
#         renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
#    output:
#         metabat2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2_bins","bin.1.fa"),
#         maxbin2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2_bins","bin.0.fa"),
#         concoct_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct_bins","bin.0.fa")
#    params:
#         sample_initial_binning_dir = os.path.join(config["output_dir"],"{sample}","initial_binning"),
#         threads = config["binning_threads"],
#         fastq_read12 = os.path.join(config["input_dir"],"{sample}*fastq")
#    conda: "utils/envs/metawrap_env.yaml"
#    shell:
#         "metawrap binning -o {params.sample_initial_binning_dir} -t {params.threads} -a {input.renamed_metagenome_assembly_file} --metabat2 --maxbin2 --concoct {params.fastq_read12}"

#rule metawrap_bin_refinement:
#    input:
#         metabat2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2_bins","bin.1.fa"),
#         maxbin2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2_bins","bin.0.fa"),
#         concoct_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct_bins","bin.0.fa")
#
#    output:
#         metabat2_bin_file = os.path.join(config["output_dir"],"{sample}","bin_refinement","metabat2_bins","bin.1.fa"),
#         maxbin2_bin_file = os.path.join(config["output_dir"],"{sample}","bin_refinement","maxbin2_bins","bin.0.fa"),
#         concoct_bin_file = os.path.join(config["output_dir"],"{sample}","bin_refinement","concoct_bins","bin.0.fa")
#    params:
#         sample_bin_refinement_dir = os.path.join(config["output_dir"],"{sample}","bin_refinement"),
#         threads = config["bin_refinement_threads"],
#         metabat2_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2_bins"),
#         maxbin2_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2_bins"),
#         concoct_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct_bins"),
#         completeness_thresh = config["completeness_thresh"],
#         contamination_thresh = config["contamination_thresh"]
#    conda: "utils/envs/metawrap_env.yaml"
#    shell:
#         "checkm data setRoot /bulk/IMCshared_bulk/shared/dbs/checkm_db;"
#         "metawrap bin_refinement -o {params.sample_bin_refinement_dir} -t {params.threads} -A {params.metabat2_bins_dir} -B {params.maxbin2_bins_dir} -C {params.concoct_bins_dir} -c {params.completeness_thresh} -x {params.contamination_thresh}"

