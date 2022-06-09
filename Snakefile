# ***************************************
# * Snakefile for metagenotate pipeline *
# ***************************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
import os

# Set the PATH environment variable for metaWRAP bin folder.
os.environ["PATH"]+=os.pathsep+"/bulk/IMCshared_bulk/shared/shared_software/metaWRAP/bin"

os.environ["GTDBTK_DATA_PATH"] = "/bulk/IMCshared_bulk/shared/dbs/gtdbtk-1.5.0/db"

SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input: 
        expand(config["output_dir"]+"/{sample}/assembly/metaspades/scaffolds.fasta",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/metaspades/long_scaffolds.fasta",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/metaspades/mapped_metaspades_assembly_reads.sam",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/metaspades/unmapped_metaspades_assembly_reads_R1.fastq",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/metaspades/unmapped_metaspades_assembly_reads_R2.fastq",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/megahit/fixed_headers_final.contigs.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/final_assembly.fasta",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/{sample}_metagenome.fasta",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/quast/transposed_report.tsv",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/prokka/{sample}_metagenome.fna",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/prokka/{sample}_metagenome.gff",sample=SAMPLES),
#        expand(config["output_dir"]+"/{sample}/assembly/metaerg/data/all.gff",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/assembly/mapped_metagenome_assembly_reads.sam",sample=SAMPLES),
   #     expand(config["output_dir"]+"/{sample}/assembly/mapped_metagenome_assembly_reads.bam",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/initial_binning/metabat2/bin.1.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/initial_binning/maxbin2_abund_list.txt",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/initial_binning/maxbin2/bin.1.fa",sample=SAMPLES),
##        expand(config["output_dir"]+"/{sample}/initial_binning/concoct/concoct_bins/bin.0.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/bin_refinement/metabat2/bin.1.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/bin_refinement/maxbin2/bin.1.fa",sample=SAMPLES),

##        expand(config["output_dir"]+"/{sample}/bin_refinement/concoct_bins/bin.0.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/bin_refinement/metawrap_" + str(config["completeness_thresh"]) + "_" + str(config["contamination_thresh"]) + "_bins.stats",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/refined_bins/{sample}_bin.1.fa",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/refined_bins/{sample}_bin.1/quast/transposed_report.tsv",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/refined_bins/{sample}_bin.1/prokka/{sample}_bin.1.fna",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/refined_bins/{sample}_bin.1/prokka/{sample}_bin.1.gff",sample=SAMPLES),
#        expand(config["output_dir"]+"/{sample}/refined_bins/{sample}_bin.1/metaerg/data/all.gff",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/refined_bins/{sample}_bin.1/checkm/checkm.tsv",sample=SAMPLES),
        expand(config["output_dir"]+"/{sample}/refined_bins/{sample}_bin.1/gtdbtk/gtdbtk.bac120.summary.tsv",sample=SAMPLES)

rule metaspades_assembly:
    input:
        fastq_read1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        fastq_read2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        metaspades_scaffolds_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/scaffolds.fasta")
    params:
        memory_in_gb = config["memory_in_gb"],
        threads = config["assembler_threads"],
        min_contig_length = config["min_contig_length"],
        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly","metaspades")
    conda: "utils/envs/metaspades_env.yaml"
    shell:
        "metaspades.py -t {params.threads} -m {params.memory_in_gb} -o {params.sample_assembly_dir} -1 {input.fastq_read1} -2 {input.fastq_read2}"

rule filter_long_scaffolds:
    input:
        metaspades_scaffolds_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/scaffolds.fasta")
    output:
        long_scaffolds_metaspades_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/long_scaffolds.fasta")
    params:
        min_sequence_length = config["min_sequence_length"]
    conda: "utils/envs/biopython_env.yaml"
    shell:
        "python utils/scripts/filter_sequences_by_length.py -i {input.metaspades_scaffolds_assembly_file} -l {params.min_sequence_length} -o {output.long_scaffolds_metaspades_file}"
        
rule map_reads_to_assembly:
    input:
        long_scaffolds_metaspades_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/long_scaffolds.fasta"),
        fastq_read1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        fastq_read2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        metaspades_scaffolds_sam_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/mapped_metaspades_assembly_reads.sam")
    params:
        threads = config["assembler_threads"],
        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly")
    conda: "utils/envs/bwa_env.yaml"
    shell:
        "bwa index {input.long_scaffolds_metaspades_file}; "
        "bwa mem -t {params.threads} {input.long_scaffolds_metaspades_file} {input.fastq_read1} {input.fastq_read2} > {output.metaspades_scaffolds_sam_file}; "

rule recover_unmapped_metaspades_reads:
    input:
        metaspades_scaffolds_sam_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/mapped_metaspades_assembly_reads.sam")
    output:
        unmapped_metaspades_scaffolds_read1_fastq_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/unmapped_metaspades_assembly_reads_R1.fastq"),
        unmapped_metaspades_scaffolds_read2_fastq_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/unmapped_metaspades_assembly_reads_R2.fastq")
    params:
        threads = config["assembler_threads"],
        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly"),
        metaspades_scaffolds_bam_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/mapped_metaspades_assembly_reads.bam"),
        unmapped_metaspades_scaffolds_bam_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/unmapped_metaspades_assembly_reads.bam")
    conda: "utils/envs/samtools_env.yaml"
    shell:
        "samtools view -S -b {input.metaspades_scaffolds_sam_file} --threads {params.threads} -o {params.metaspades_scaffolds_bam_file}; "
        "samtools view -b -f 4 {params.metaspades_scaffolds_bam_file} --threads {params.threads} -o {params.unmapped_metaspades_scaffolds_bam_file}; "
        "samtools fastq {params.unmapped_metaspades_scaffolds_bam_file} --threads {params.threads} -1 {output.unmapped_metaspades_scaffolds_read1_fastq_file} -2 {output.unmapped_metaspades_scaffolds_read2_fastq_file}; "

rule megahit_unmapped_metaspades_reads:
    input:
        unmapped_metaspades_scaffolds_read1_fastq_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/unmapped_metaspades_assembly_reads_R1.fastq"),
        unmapped_metaspades_scaffolds_read2_fastq_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/unmapped_metaspades_assembly_reads_R2.fastq")
    output:
        fixed_megahit_final_contigs_file = os.path.join(config["output_dir"],"{sample}","assembly","megahit/fixed_headers_final.contigs.fa")
    params:
        megahit_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","megahit/final.contigs.fa"),
        threads = config["assembler_threads"],
        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly"),
        memory_in_gb = config["memory_in_gb"]
    conda: "utils/envs/megahit_env.yaml"
    shell:
        "megahit -1 {input.unmapped_metaspades_scaffolds_read1_fastq_file} -2 {input.unmapped_metaspades_scaffolds_read2_fastq_file} -t {params.threads} -m {params.memory_in_gb}000000000 -f -o {params.sample_assembly_dir}/megahit; "
        "sed 's/ /_/g' {params.megahit_assembly_file} > {output.fixed_megahit_final_contigs_file} "

rule combine_and_sort_assembly:
    input:
        long_scaffolds_metaspades_file = os.path.join(config["output_dir"],"{sample}","assembly","metaspades/long_scaffolds.fasta"),
        fixed_megahit_final_contigs_file = os.path.join(config["output_dir"],"{sample}","assembly","megahit/fixed_headers_final.contigs.fa")
    output:
        metagenome_final_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","final_assembly.fasta")
    params:
        threads = config["assembler_threads"],
        combined_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","combined_file.fasta"),
        long_scaffolds_megahit_file = os.path.join(config["output_dir"],"{sample}","assembly","megahit/fixed_headers_long.contigs.fa"),
        min_sequence_length = config["min_sequence_length"],
        memory_in_gb = config["memory_in_gb"]

    conda: "utils/envs/biopython_env.yaml"
    shell:
        "python utils/scripts/filter_sequences_by_length.py -i {input.fixed_megahit_final_contigs_file} -l {params.min_sequence_length} -o {params.long_scaffolds_megahit_file}; "
        "cat {input.long_scaffolds_metaspades_file} {params.long_scaffolds_megahit_file} > {params.combined_assembly_file}; "
        "python utils/scripts/sort_assembly_by_length.py -i {params.combined_assembly_file} -o {output.metagenome_final_assembly_file}; "
        
rule rename_final_assembly_file:
    input:
        metagenome_final_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","final_assembly.fasta")
    output:
        renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
    conda: "utils/envs/biopython_env.yaml"
    shell:
        "python utils/scripts/fix_fasta_header_length.py -i {input.metagenome_final_assembly_file} -o {output.renamed_metagenome_assembly_file}"

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
        prokka_fna_file = os.path.join(config["output_dir"],"{sample}","assembly","prokka","{sample}_metagenome.fna"),
        prokka_gff_file = os.path.join(config["output_dir"],"{sample}","assembly","prokka","{sample}_metagenome.gff")
    params:
        prokka_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly","prokka"),
        threads = config["prokka_threads"],
	prefix = "{sample}_metagenome"
    conda: "utils/envs/prokka_env.yaml"
    shell:
       "prokka --metagenome --outdir {params.prokka_assembly_dir} --prefix {params.prefix} {input.renamed_metagenome_assembly_file} --cpus {params.threads} --rfam 1 --force"


#rule metaerg_assembly:
#    input:
#        renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
#    output:
##        metaerg_fna_file = os.path.join(config["output_dir"],"{sample}","assembly","metaerg","{sample}_metagenome.fna"),
#        metaerg_gff_file = os.path.join(config["output_dir"],"{sample}","assembly","metaerg","data","all.gff")
#    params:
#        assembly_metaerg_dir = os.path.join(config["output_dir"],"{sample}","assembly","metaerg"),
#        metaerg_database_path = config["metaerg_database_path"],
#        locustag = "{sample}_metagenome",
#        threads = config["metaerg_threads"]
#    shell:
#       "singularity run -H $HOME -B {params.metaerg_database_path}:/NGStools/metaerg/db -B /work:/work -B /bulk:/bulk /global/software/singularity/images/software/metaerg2.sif /NGStools/metaerg/bin/metaerg.pl --mincontiglen 200 --gcode 11 --gtype meta --minorflen 180 --cpus {params.threads} --evalue 1e-05 --identity 20 --coverage 70 --locustag {params.locustag} --force --outdir {params.assembly_metaerg_dir} {input.renamed_metagenome_assembly_file}"

rule map_reads_to_metagenome:
    input:
        renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta"),
        fastq_read1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        fastq_read2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        metagenome_sam_file = os.path.join(config["output_dir"],"{sample}","assembly","mapped_metagenome_assembly_reads.sam")
    params:
        threads = config["binning_threads"],
        sample_assembly_dir = os.path.join(config["output_dir"],"{sample}","assembly")
    conda: "utils/envs/bwa_env.yaml"
    shell:
        "bwa index {input.renamed_metagenome_assembly_file}; "
        "bwa mem -t {params.threads} {input.renamed_metagenome_assembly_file} {input.fastq_read1} {input.fastq_read2} > {output.metagenome_sam_file}; "


rule sort_mapped_metagenome_reads:
    input:
        metagenome_sam_file = os.path.join(config["output_dir"],"{sample}","assembly","mapped_metagenome_assembly_reads.sam")
    output:
        metagenome_bam_file = os.path.join(config["output_dir"],"{sample}","assembly","mapped_metagenome_assembly_sorted_reads.bam")
    params:
        threads = config["binning_threads"],
    conda: "utils/envs/samtools_env.yaml"
    shell:
        "samtools sort -@ {params.threads} -O BAM -o {output.metagenome_bam_file} {input.metagenome_sam_file}; "

rule metabat2_binning:
    input:
         metagenome_bam_file = os.path.join(config["output_dir"],"{sample}","assembly","mapped_metagenome_assembly_sorted_reads.bam"),
         renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
    output:
         metabat2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2","bin.1.fa"),
         maxbin2_abund_list_file = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2_abund_list.txt"),
         metabat2_checkpoint_file = os.path.join(config["output_dir"],"{sample}","initial_binning","working_dir","metabat2","metabat2_checkpoint.txt"),
    params:
         metabat2_depth_file = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2_depth.txt"),
         maxbin2_depth_file = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2_depth.txt"),
         maxbin2_abund_file = os.path.join(config["output_dir"],"{sample}","initial_binning","{sample}_maxbin2_abund.txt"),
         metabat2_bin_outfile_prefix = os.path.join(config["output_dir"],"{sample}","initial_binning","working_dir","metabat2","bin"),
         metabat2_working_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","working_dir","metabat2"),
         metabat2_bin_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2"),         
         threads = config["binning_threads"],
         min_sequence_length = config["min_sequence_length"],
    conda: "utils/envs/metabat2_env.yaml"
    shell:
         "jgi_summarize_bam_contig_depths --outputDepth {params.metabat2_depth_file} {input.metagenome_bam_file}; "
         "mkdir -p {params.metabat2_working_dir}; "
         "metabat2 -i {input.renamed_metagenome_assembly_file} -a {params.metabat2_depth_file} -o {params.metabat2_bin_outfile_prefix} -m {params.min_sequence_length} -t {params.threads} --unbinned; "
         "jgi_summarize_bam_contig_depths --outputDepth {params.maxbin2_depth_file} --noIntraDepthVariance {input.metagenome_bam_file}; "
         "tail -n+2 {params.maxbin2_depth_file} | cut -f1,3 > {params.maxbin2_abund_file}; "
         "echo \"{params.maxbin2_abund_file}\" > {output.maxbin2_abund_list_file}; "
         "mkdir -p {params.metabat2_bin_dir}; "
         "cp {params.metabat2_bin_outfile_prefix}.[0-9]*.fa {params.metabat2_bin_dir}; "
         "echo \"metabat2_binning rule completed. Done!\" > {output.metabat2_checkpoint_file}; "

rule maxbin2_binning:
    input:
         metagenome_bam_file = os.path.join(config["output_dir"],"{sample}","assembly","mapped_metagenome_assembly_sorted_reads.bam"),
         renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta"),
         maxbin2_abund_list_file = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2_abund_list.txt"),
    output:
         maxbin2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2","bin.1.fa"),
         maxbin2_checkpoint_file = os.path.join(config["output_dir"],"{sample}","initial_binning","working_dir","maxbin2","maxbin2_checkpoint.txt"),
    params:
         maxbin2_bin_outfile_prefix = os.path.join(config["output_dir"],"{sample}","initial_binning","working_dir","maxbin2","bin"),
         maxbin2_working_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","working_dir","maxbin2"),
         maxbin2_bin_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2"),
         maxbin2_path = config["maxbin2_path"], 
         threads = config["binning_threads"],
         min_sequence_length = config["min_sequence_length"],
         markers = config["maxbin2_markers"],
    conda: "utils/envs/maxbin2_env.yaml"
    shell:
         "mkdir -p {params.maxbin2_working_dir}; "         
         "perl {params.maxbin2_path}/run_MaxBin.pl -contig {input.renamed_metagenome_assembly_file} -markerset {params.markers} -thread {params.threads} -min_contig_length {params.min_sequence_length} -out {params.maxbin2_bin_outfile_prefix} -abund_list {input.maxbin2_abund_list_file}; "
         "mkdir -p {params.maxbin2_bin_dir}; "
         "for bin_file in $(ls {params.maxbin2_working_dir} | grep \"\.fasta\"); "
         "do echo $bin_file; filename=$(basename $bin_file '.fasta'); bin_count=$(echo $filename | sed 's/bin\.0\+//g'); echo $bin_count; new_filename=\"bin.$bin_count.fa\"; echo $new_filename; cp {params.maxbin2_working_dir}/$bin_file {params.maxbin2_bin_dir}/$new_filename; done; "
         "echo \"maxbin2_binning rule completed. Done!\" > {output.maxbin2_checkpoint_file}; "
 
#rule metawrap_concoct_binning:
#    input:
#         renamed_metagenome_assembly_file = os.path.join(config["output_dir"],"{sample}","assembly","{sample}_metagenome.fasta")
#    output:
#         concoct_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct","concoct_bins","bin.0.fa")
#    params:
#         metawrap_path = config["metawrap_path"],
#         sample_initial_binning_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct"),
#         threads = config["binning_threads"],
#         fastq_read12 = os.path.join(config["input_dir"],"{sample}*fastq")
#    conda: "utils/envs/metawrap_env.yaml"
#    shell:
#         "{params.metawrap_path}/metawrap binning -o {params.sample_initial_binning_dir} -t {params.threads} -a {input.renamed_metagenome_assembly_file} --concoct {params.fastq_read12}"

rule metawrap_bin_refinement:
    input:
         metabat2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2","bin.1.fa"),
         maxbin2_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2","bin.1.fa"),
##         concoct_bin_file = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct","concoct_bins","bin.0.fa")

    output:
         metabat2_bin_file = os.path.join(config["output_dir"],"{sample}","bin_refinement","metabat2_bins","bin.1.fa"),
         maxbin2_bin_file = os.path.join(config["output_dir"],"{sample}","bin_refinement","maxbin2_bins","bin.1.fa"),
##         concoct_bin_file = os.path.join(config["output_dir"],"{sample}","bin_refinement","concoct_bins","bin.0.fa"),
         metawrap_refine_bin_stats = os.path.join(config["output_dir"],"{sample}","bin_refinement","_".join(["metawrap",str(config["completeness_thresh"]),str(config["contamination_thresh"]),"bins.stats"])),
         refined_bin_file = os.path.join(config["output_dir"],"{sample}","bin_refinement","_".join(["metawrap",str(config["completeness_thresh"]),str(config["contamination_thresh"]),"bins"]), "bin.1.fa"),
    params:
         metawrap_path = config["metawrap_path"],
         checkm_database = config["checkm_database_path"],
         sample_bin_refinement_dir = os.path.join(config["output_dir"],"{sample}","bin_refinement"),
         threads = config["bin_refinement_threads"],
         metabat2_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","metabat2"),
         maxbin2_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","maxbin2"),
##         concoct_bins_dir = os.path.join(config["output_dir"],"{sample}","initial_binning","concoct","concoct_bins"),
         completeness_thresh = config["completeness_thresh"],
         contamination_thresh = config["contamination_thresh"]
    conda: "utils/envs/metawrap_bin_refinement_env.yaml"
    shell:
         "echo \"{params.checkm_database}\" | checkm data setRoot; "
##         "{params.metawrap_path}/metawrap bin_refinement -o {params.sample_bin_refinement_dir} -t {params.threads} -A {params.metabat2_bins_dir} -B {params.maxbin2_bins_dir} -C {params.concoct_bins_dir} -c {params.completeness_thresh} -x {params.contamination_thresh}"
         "{params.metawrap_path}/metawrap bin_refinement -o {params.sample_bin_refinement_dir} -t {params.threads} -A {params.metabat2_bins_dir} -B {params.maxbin2_bins_dir} -c {params.completeness_thresh} -x {params.contamination_thresh}"

## Going to recreate this so just a directory is used. I can usethe find command.
rule rename_refined_bin_file:
    input:
        refined_bin_file = os.path.join(config["output_dir"],"{sample}","bin_refinement","_".join(["metawrap",str(config["completeness_thresh"]),str(config["contamination_thresh"]),"bins"]), "bin.1.fa")
    output:
        renamed_refined_bin_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1.fa")
    params:
        metawrap_bin_refinement_dir = os.path.join(config["output_dir"],"{sample}","bin_refinement","_".join(["metawrap",str(config["completeness_thresh"]),str(config["contamination_thresh"]),"bins"])),
        refined_bins_dir = os.path.join(config["output_dir"],"{sample}","refined_bins"),
        sample_name = "{sample}",
    shell:
       "COUNTER=1; "
       "for bin_file in $(ls {params.metawrap_bin_refinement_dir} | grep \"\.fa\"); "
       "do echo $bin_file; "
       "filename=$(basename $bin_file \".fa\"); "
       "renamed_refined_bin_file=\"{params.refined_bins_dir}/{params.sample_name}_bin.$COUNTER.fa\"; "
       "cp {params.metawrap_bin_refinement_dir}/$bin_file $renamed_refined_bin_file; "
       "echo \"$bin_file\t{params.sample_name}_bin.$COUNTER.fa\" >> {params.refined_bins_dir}/{params.sample_name}_bin_rename_index.txt; "
       "COUNTER=$((COUNTER+1)); "
       "done"

rule quast_refined_bins:
     input:
        refined_bin_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1.fa")
     output:
        quast_transposed_report_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1", "quast","transposed_report.tsv")
     params:
        refined_bins_dir = os.path.join(config["output_dir"],"{sample}","refined_bins"),
        threads = config["quast_threads"]
     conda: "utils/envs/quast_env.yaml"
     shell:
       "for bin_file in $(ls {params.refined_bins_dir} | grep \"\.fa\"); "
       "do echo $bin_file; "
       "filename=$(basename $bin_file \".fa\"); "
       "bin_dir=\"{params.refined_bins_dir}/$filename\"; "
       "mkdir -p $bin_dir; "
       "quast_bin_dir=\"$bin_dir/quast\"; "
       "mkdir -p $quast_bin_dir; "
       "refined_bin_file=\"{params.refined_bins_dir}/$bin_file\"; "
       "quast.py --output-dir $quast_bin_dir --threads {params.threads} $refined_bin_file; "
       "done"

rule prokka_refined_bins:
    input:
        refined_bin_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1.fa")
    output:
        prokka_fna_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1","prokka","{sample}_bin.1.fna"),
        prokka_gff_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1","prokka","{sample}_bin.1.gff")
    params:
        refined_bins_dir = os.path.join(config["output_dir"],"{sample}","refined_bins"),
        threads = config["prokka_threads"],
    conda: "utils/envs/prokka_env.yaml"
    shell:
       "for bin_file in $(ls {params.refined_bins_dir} | grep \"\.fa\"); "
       "do echo $bin_file; "
       "filename=$(basename $bin_file \".fa\"); "
       "bin_dir=\"{params.refined_bins_dir}/$filename\"; "
       "mkdir -p $bin_dir; "
       "prokka_bin_dir=\"$bin_dir/prokka\"; "
       "mkdir -p $prokka_bin_dir; "
       "refined_bin_file=\"{params.refined_bins_dir}/$bin_file\"; "
       "prefix=$filename; "
       "prokka --metagenome --outdir $prokka_bin_dir --prefix $prefix $refined_bin_file --cpus {params.threads} --rfam 1 --force; "
       "done"

#rule metaerg_refined_bins:
#    input:
#        refined_bin_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1.fa")
#    output:
###        metaerg_fna_file = os.path.join(config["output_dir"],"{sample}","assembly","metaerg","{sample}_metagenome.fna"),
#        metaerg_gff_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1","metaerg","data","all.gff")
#    params:
#        metaerg_database_path = config["metaerg_database_path"],
#        refined_bins_dir = os.path.join(config["output_dir"],"{sample}","refined_bins"),
#        threads = config["metaerg_threads"]
#    shell:
#       "for bin_file in $(ls {params.refined_bins_dir} | grep \"\.fa\"); "
#       "do echo $bin_file; "
#       "filename=$(basename $bin_file \".fa\"); "
#       "bin_dir=\"{params.refined_bins_dir}/$filename\"; "
#       "mkdir -p $bin_dir; "
#       "metaerg_bin_dir=\"$bin_dir/metaerg\"; "
#       "mkdir -p $metaerg_bin_dir; "
#       "refined_bin_file=\"{params.refined_bins_dir}/$bin_file\"; "
#       "locus_tag=$filename; "
#       "singularity run -H $HOME -B {params.metaerg_database_path}:/NGStools/metaerg/db -B /work:/work -B /bulk:/bulk /global/software/singularity/images/software/metaerg2.sif /NGStools/metaerg/bin/metaerg.pl --mincontiglen 200 --gcode 11 --gtype meta --minorflen 180 --cpus {params.threads} --evalue 1e-05 --identity 20 --coverage 70 --locustag $locus_tag --force --outdir $metaerg_bin_dir $refined_bin_file; "
#       "done"

rule checkm_refined_bins:
    input:
        refined_bin_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1.fa")
    output:
        checkm_table_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1","checkm","checkm.tsv")
    params:
        checkm_database = config["checkm_database_path"],
        refined_bins_dir = os.path.join(config["output_dir"],"{sample}","refined_bins"),
        threads = config["checkm_threads"]
    conda: "utils/envs/checkm_env.yaml"
    shell:
       "echo \"{params.checkm_database}\" | checkm data setRoot; "
       "for bin_file in $(ls {params.refined_bins_dir} | grep \"\.fa\"); "
       "do echo $bin_file; "
       "filename=$(basename $bin_file \".fa\"); "
       "bin_dir=\"{params.refined_bins_dir}/$filename\"; "
       "mkdir -p $bin_dir; "
       "checkm_bin_dir=\"$bin_dir/checkm\"; "
       "mkdir -p $checkm_bin_dir; "
       "checkm_table_file=\"$checkm_bin_dir/checkm.tsv\"; "
       "cp {params.refined_bins_dir}/$bin_file $checkm_bin_dir/$bin_file; "
       "checkm lineage_wf -t {params.threads} -x fa --tab_table --file $checkm_table_file $checkm_bin_dir $checkm_bin_dir; "
       "done"

rule gtdbtk_refined_bins:
    input:
        refined_bin_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1.fa")
    output:
        gtdbtk_refined_bin_file = os.path.join(config["output_dir"],"{sample}","refined_bins","{sample}_bin.1","gtdbtk","gtdbtk.bac120.summary.tsv")
    params:
       refined_bins_dir = os.path.join(config["output_dir"],"{sample}","refined_bins"),
       gtdbtk_data_path = config["gtdbtk_database_path"],
       threads = config["gtdbtk_threads"]
    conda: "utils/envs/gtdbtk_env.yaml"
    shell:
       "GTDBTK_DATA_PATH=\"{params.gtdbtk_data_path}\"; "
       "for bin_file in $(ls {params.refined_bins_dir} | grep \"\.fa\"); "
       "do echo $bin_file; "
       "filename=$(basename $bin_file \".fa\"); "
       "bin_dir=\"{params.refined_bins_dir}/$filename\"; "
       "mkdir -p $bin_dir; "
       "gtdbtk_bin_dir=\"$bin_dir/gtdbtk\"; "
       "mkdir -p $gtdbtk_bin_dir; "
       "cp {params.refined_bins_dir}/$bin_file $gtdbtk_bin_dir/$bin_file; "
       "gtdbtk classify_wf --genome_dir $gtdbtk_bin_dir --extension \"fa\" --cpus {params.threads} --out_dir $gtdbtk_bin_dir; "
       "done"

