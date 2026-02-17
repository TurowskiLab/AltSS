import os, sys, psutil
import pandas as pd

########## CONFIGURATION ##########
# Check for required config values
required_params = ["raw_path", "gtf", "fasta", "outdir", "read_length"]
missing_params = [p for p in required_params if not config.get(p)]

if missing_params:
    print(f"ERROR: Missing required parameters: {', '.join(missing_params)}")
    print("Please provide them via config.yml file")
    sys.exit(1)

RAW_PATH = config["raw_path"]
GTF = config["gtf"]
FASTA = config["fasta"]
OUTDIR = config["outdir"]
if not OUTDIR.endswith("/"):
    OUTDIR += "/"
if not RAW_PATH.endswith("/"):
    RAW_PATH += "/"
READLENGTH = config["read_length"]
READTYPE = config.get("read_type", "paired")
CONTROL = config.get("control_name", "control")
NTHREAD = config.get("nthread", 8)
OVERHANG = READLENGTH - 1
TOOL_TYPE = "rmats"

if READTYPE == "paired":
    FASTQ_PATTERN = "_1.fastq.gz"
else:
    FASTQ_PATTERN = ".fastq.gz"

print(f"Using {NTHREAD} threads")
##################################

# Ensure output directory exists
os.makedirs(OUTDIR, exist_ok=True)

########## PARSE SAMPLES FROM FILES ##########

fastq_files = [f for f in os.listdir(RAW_PATH) if f.endswith(FASTQ_PATTERN)]

if not fastq_files:
    raise ValueError(f"No FASTQ files (*{FASTQ_PATTERN}) found in {RAW_PATH}")

longNames = [f.replace(FASTQ_PATTERN, "") for f in fastq_files]

df_names = pd.DataFrame({
    'longName': longNames,
    'sample': ["_".join(n.split("_")[:3]) for n in longNames],
    'bam_name': ["_".join(n.split("_")[:-1]) for n in longNames],
    'rmats_run': ["_".join(n.split("_")[:2]) for n in longNames],
    'control': [n.split("_")[0] + "_" + CONTROL for n in longNames]
})

# Export metadata for reference
df_names.to_csv(OUTDIR + 'metadata.tsv', sep="\t", index=False)
print(df_names)

# Extract variables
SAMPLES = df_names['sample'].unique().tolist()
bamlist = df_names['bam_name'].unique().tolist()
ctrl_dict = df_names[['rmats_run', 'control']].drop_duplicates().set_index('rmats_run')['control'].to_dict()
rmats_run = df_names[~df_names['rmats_run'].str.contains(CONTROL)]['rmats_run'].unique().tolist()
print(bamlist)
print(rmats_run)
########## OUTPUTS ##########

rule all:
        input:
                OUTDIR + "index.done",
                expand(OUTDIR + "alignment/{sample}.bam",sample=SAMPLES),
                expand(OUTDIR + "alignment/{sample}.bam.bai",sample=SAMPLES),
                OUTDIR + "cleaninig.done",
                OUTDIR + "FeatureCounts/featureCounts_multimappers.list",
                OUTDIR + "FeatureCounts/featureCounts_uniq.list",
                expand(OUTDIR + "BigWig/{sample}_raw_plus.bw",sample=SAMPLES),
                expand(OUTDIR + "BigWig/{sample}_raw_minus.bw",sample=SAMPLES),
                expand(OUTDIR + "BigWig/{sample}_CPM_plus.bw",sample=SAMPLES),
                expand(OUTDIR + "BigWig/{sample}_CPM_minus.bw",sample=SAMPLES),
                expand(OUTDIR + "bam_lists/bam_{bamlist}.txt", bamlist=bamlist),
                expand(OUTDIR + "rmats_{rmats_run}.done", rmats_run=rmats_run)

########## ALIGNMENT ##########

rule generate_star_index:
        input:
                fasta = FASTA,
                gtf = GTF
        output:
                touch(OUTDIR + 'index.done')
        params:
                index_dir = OUTDIR + "STAR_index/",
                overhang = OVERHANG,
                threads = NTHREAD
        conda:
                "processingPipeline"
        shell:
                """
                mkdir -p {params.index_dir} &&
                STAR --runMode genomeGenerate \
                --genomeDir {params.index_dir} \
                --genomeFastaFiles {input.fasta} \
                --sjdbGTFfile {input.gtf} \
                --sjdbOverhang {params.overhang} \
                --runThreadN {params.threads}
                """

rule load_genome:
        input:
                index_check = OUTDIR + 'index.done',
        params:
                index_dir = OUTDIR + "STAR_index/"
        output:
                touch(OUTDIR + 'loading.done')
        conda:
                "processingPipeline"
        shell:
                "STAR --genomeLoad LoadAndExit --genomeDir {params.index_dir}"


def get_fastq_inputs(wildcards):
    """Return appropriate FASTQ inputs based on read type"""
    if READTYPE == "paired":
        return {
            "read1": RAW_PATH + f"{wildcards.sample}_1.fastq.gz",
            "read2": RAW_PATH + f"{wildcards.sample}_2.fastq.gz"
        }
    else:  # single
        return {
            "read1": RAW_PATH + f"{wildcards.sample}.fastq.gz"
        }

rule align:
        input:
                unpack(get_fastq_inputs),
                idx = OUTDIR + 'loading.done'
        params:
                index_dir = OUTDIR + "STAR_index/",
                prefix = OUTDIR + "alignment/{sample}_STAR_",
                reads = lambda wildcards, input: f"{input.read1} {input.read2}" if READTYPE == "paired" else input.read1
        output:
                bam = OUTDIR + "alignment/{sample}_STAR_Aligned.out.bam"
        conda:
                "processingPipeline"
        shell:
                "STAR --outFileNamePrefix {params.prefix} --readFilesCommand zcat --genomeDir {params.index_dir} --genomeLoad LoadAndKeep --outSAMtype BAM Unsorted --readFilesIn {params.reads}"

########## POSTPROCESSING ##########

rule sort:
        input:
                bam = OUTDIR + "alignment/{sample}_STAR_Aligned.out.bam"
        output:
                bam = OUTDIR + "alignment/{sample}.bam",
                bai = OUTDIR + "alignment/{sample}.bam.bai"
        conda:
                "processingPipeline"
        shell:
                """
                samtools sort {input.bam} > {output.bam}
                samtools index {output.bam}
                """

rule clean:
      input:
              expand(OUTDIR + "alignment/{sample}.bam.bai",sample=SAMPLES)
      output:
              touch(OUTDIR + "cleaninig.done")
      shell:
              """
              rm Aligned.out.sam Log.final.out Log.out Log.progress.out SJ.out.tab
              """

rule featureCounts:
        input:
                bam = expand(OUTDIR + "alignment/{sample}.bam",sample=SAMPLES) #use list of files
        output:
                multi = OUTDIR + "FeatureCounts/featureCounts_multimappers.list",
                uniq = OUTDIR + "FeatureCounts/featureCounts_uniq.list"
        params:
                gtf=GTF
        conda:
                "processingPipeline"
        shell:
                """
                featureCounts -M -s 1 -a {params.gtf} -o {output.multi} {input.bam}
                featureCounts -s 1 -a {params.gtf} -o {output.uniq} {input.bam}
                """

rule BigWigs_CPM:
        input:
                bam = OUTDIR + "alignment/{sample}.bam",
                bai = OUTDIR + "alignment/{sample}.bam.bai"
        output:
                bwP = OUTDIR + "BigWig/{sample}_CPM_plus.bw",
                bwM = OUTDIR + "BigWig/{sample}_CPM_minus.bw"
        conda:
                "processingPipeline"
        shell:
                """
                bamCoverage --bam {input.bam} -of bigwig -o {output.bwP} --filterRNAstrand reverse --normalizeUsing CPM --binSize 1
                bamCoverage --bam {input.bam} -of bigwig -o {output.bwM} --filterRNAstrand forward --normalizeUsing CPM --binSize 1
                """

rule BigWigs_raw:
        input:
                bam = OUTDIR + "alignment/{sample}.bam",
                bai = OUTDIR + "alignment/{sample}.bam.bai"
        output:
                bwP = OUTDIR + "BigWig/{sample}_raw_plus.bw",
                bwM = OUTDIR + "BigWig/{sample}_raw_minus.bw"
        conda:
                "processingPipeline"
        shell:
                """
                bamCoverage --bam {input.bam} -of bigwig -o {output.bwP} --filterRNAstrand reverse --binSize 1
                bamCoverage --bam {input.bam} -of bigwig -o {output.bwM} --filterRNAstrand forward --binSize 1
                """

########## RMATS ##########

rule list_bams:
        input:
                expand(OUTDIR + "alignment/{sample}.bam",sample=SAMPLES),
                expand(OUTDIR + "alignment/{sample}.bam.bai",sample=SAMPLES),
        output:
                bl = OUTDIR + "bam_lists/bam_{bamlist}.txt"
        params:
                bn = "{bamlist}",
                alignment_dir = OUTDIR + "alignment",
                bam_lists_dir = OUTDIR + "bam_lists",
                tool_type = TOOL_TYPE
        shell:
                "scripts/list_bams.sh {params.bn} {params.alignment_dir} {params.bam_lists_dir} {params.tool_type}"

def get_ctrl(wildcards):
    rmats_run = wildcards.rmats_run
    return OUTDIR+"bam_lists/bam_"+ctrl_dict[rmats_run]+".txt"


rule run_rmats:
        input:
                b2=OUTDIR+"bam_lists/bam_{rmats_run}.txt",
                b1=get_ctrl,
        output:
                touch(OUTDIR+"rmats_{rmats_run}.done")
        params:
                readLength=READLENGTH,
                nthread=NTHREAD,
                od=OUTDIR + "rmats_{rmats_run}/",
                gtf=GTF,
                readType=READTYPE
        conda:
                "rmats"
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
                --od {params.od} \
                --tmp {params.od} \
                -t {params.readType}
                """