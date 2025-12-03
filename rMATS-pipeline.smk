import re, os, subprocess
import pandas as pd

########## CHANGE THESE ##########
raw_path = "/home/rgoitia/00_raw_data/01_Choquet_POLR3-HLD/" #Path to raw RNA-seq data
#references
STAR_INDEX = "/home/rgoitia/01_reference_sequences/hg41/hg41_STAR_index/" #Path to STAR index
GTF = "/home/rgoitia/01_reference_sequences/hg41/hg41_annotation_gencode_tRNA_rRNA.gtf" #Path to reference GTF file
##################################

#parsing file names and preparatory jobs
name_elem = "_1.fastq.gz"  #write file ending here
longName = [n.replace(name_elem,"") for n in os.listdir(raw_path) if n.endswith(name_elem)]
SAMPLES = ["_".join(n.split("_")[:3]) for n in longName]

df_names = pd.DataFrame({
        'longName' : longName,
        'bam_name': ["_".join(n.split("_")[:-1]) for n in longName],
        'control'         : [n.split("_")[0] + "_wt" for n in longName],
        'rmats_run'       : ["_".join(n.split("_")[:2]) for n in longName],
        'name'    : SAMPLES}).set_index("name")

print(df_names)
df_names.to_csv('names.tab', sep="\t")

bamlist = df_names['bam_name'].unique()
ctrl_dict = df_names[['control','rmats_run']].drop_duplicates().set_index('rmats_run')['control'].to_dict()
rmats_run = df_names.loc[~df_names['rmats_run'].str.contains('wt'), 'rmats_run'].unique()

print(bamlist)
print(ctrl_dict)
print(rmats_run)

#SnakeMake pipeline

########## OUTPUTS ##########

rule all:
        input:
                expand("alignment/{sample}.bam",sample=SAMPLES),
                expand("alignment/{sample}.bam.bai",sample=SAMPLES),
                # "cleaninig.done",
                "FeatureCounts/featureCounts_multimappers.list",
                "FeatureCounts/featureCounts_uniq.list",
                expand("BigWig/{sample}_raw_plus.bw",sample=SAMPLES),
                expand("BigWig/{sample}_raw_minus.bw",sample=SAMPLES),
                expand("BigWig/{sample}_CPM_plus.bw",sample=SAMPLES),
                expand("BigWig/{sample}_CPM_minus.bw",sample=SAMPLES),
                expand("bam_lists/bam_{bamlist}.txt", bamlist=bamlist),
                expand("rmats_output/{rmats_run}/",rmats_run=rmats_run)


########## PREPROCESSING ##########

########## ALIGNMENT ##########

rule load_genome:
        input:
                index_check = STAR_INDEX + "SAindex",
        params:
                index_dir = STAR_INDEX
        output:
                touch('loading.done')
        conda:
                "envs/processing.yml"
        shell:
                "STAR --genomeLoad LoadAndExit --genomeDir {params.index_dir}"

rule align:
        input:
                read1 = raw_path + "{sample}_1.fastq.gz",
                read2 = raw_path + "{sample}_2.fastq.gz",
                idx = 'loading.done'
        params:
                index_dir = STAR_INDEX,
                prefix = "alignment/{sample}_STAR_"
        output:
                bam = "alignment/{sample}_STAR_Aligned.out.bam"
        conda:
                "envs/processing.yml"
        shell:
                "STAR --outFileNamePrefix {params.prefix} --readFilesCommand zcat --genomeDir {params.index_dir} --genomeLoad LoadAndKeep --outSAMtype BAM Unsorted --readFilesIn {input.read1} {input.read2}"

########## POSTPROCESSING ##########

rule sort:
        input:
                bam = "alignment/{sample}_STAR_Aligned.out.bam"
        output:
                bam = "alignment/{sample}.bam",
                bai = "alignment/{sample}.bam.bai"
        conda:
                "envs/processing.yml"
        shell:
                """
                samtools sort {input.bam} > {output.bam}
                samtools index {output.bam}
                """

# rule clean:
#       input:
#               expand("alignment/{sample}.bam.bai",sample=SAMPLES)
#       output:
#               touch("cleaninig.done")
#       conda:
#               "envs/processing.yml"
#       shell:
#               """
#               rm -r logs/
#               rm Aligned.out.sam Log.final.out Log.out Log.progress.out SJ.out.tab
#               """

rule featureCounts:
        input:
                bam = expand("alignment/{sample}.bam",sample=SAMPLES) #use list of files
        output:
                multi = "FeatureCounts/featureCounts_multimappers.list",
                uniq = "FeatureCounts/featureCounts_uniq.list"
        params:
                gtf=GTF
        conda:
                "envs/processing.yml"
        shell:
                """
                featureCounts -M -s 1 -a {params.gtf} -o {output.multi} {input.bam}
                featureCounts -s 1 -a {params.gtf} -o {output.uniq} {input.bam}
                """

rule BigWigs_CPM:
        input:
                bam = "alignment/{sample}.bam",
                bai = "alignment/{sample}.bam.bai"
        output:
                bwP = "BigWig/{sample}_CPM_plus.bw",
                bwM = "BigWig/{sample}_CPM_minus.bw"
        conda:
                "envs/processing.yml"
        shell:
                """
                bamCoverage --bam {input.bam} -of bigwig -o {output.bwP} --filterRNAstrand reverse --normalizeUsing CPM --binSize 1
                bamCoverage --bam {input.bam} -of bigwig -o {output.bwM} --filterRNAstrand forward --normalizeUsing CPM --binSize 1
                """

rule BigWigs_raw:
        input:
                bam = "alignment/{sample}.bam",
                bai = "alignment/{sample}.bam.bai"
        output:
                bwP = "BigWig/{sample}_raw_plus.bw",
                bwM = "BigWig/{sample}_raw_minus.bw"
        conda:
                "envs/processing.yml"
        shell:
                """
                bamCoverage --bam {input.bam} -of bigwig -o {output.bwP} --filterRNAstrand reverse --binSize 1
                bamCoverage --bam {input.bam} -of bigwig -o {output.bwM} --filterRNAstrand forward --binSize 1
                """

########## RMATS ##########

rule list_bams:
        input:
                expand("alignment/{sample}.bam",sample=SAMPLES),
                expand("alignment/{sample}.bam.bai",sample=SAMPLES),
        output:
                bl = "bam_lists/bam_{bamlist}.txt"
        params:
                bn = "{bamlist}"
        shell:
                "scripts/list_bams.sh {params.bn} alignment bam_lists"

def get_ctrl(wildcards):
    rmats_run = wildcards.rmats_run
    return "bam_lists/bam_"+ctrl_dict[rmats_run]+".txt"

rule run_rmats:
        input:
                b2="bam_lists/bam_{rmats_run}.txt",
                b1=get_ctrl,
        output:
                directory("rmats_output/{rmats_run}/")
        params:
                readLength=50,
                nthread=4,
                od="rmats_output/",
                gtf=GTF
        conda:
                "envs/rmats.yml"
        shell:
                """
                python \
                rmats-turbo/rmats.py \
                --b1 {input.b1} \
                --b2 {input.b2} \
                --gtf {params.gtf} \
                --readLength {params.readLength} \
                --variable-read-length \
                --nthread {params.nthread} \
                --od {output} \
                --tmp {output}
                """
