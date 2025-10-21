# rMATS-pipeline

Pipeline from RNA-Seq data to rMATS alternative splicing events

## Introduction

This repository contains instructions to run rMATS Turbo (refer to <https://github.com/Xinglab/rmats-turbo>) on RNA-seq data after alignment. This is currently utilizing BAM files from different sample groups and a GTF annotation file to identify alternative splicing events.

## Requirements

Tested in WSL2 and linux server

- Conda 24.11.3

## 1. Prepare rMATS environment

1. Download `rmats-turbo` in a folder with the same name.

    ```bash
    git clone https://github.com/Xinglab/rmats-turbo.git
    ```

2. Create conda environment for snakemake

    ```bash
    conda env create -n snakemake -f envs/snakemake.yml
    ```

3. Activate rmats environment, and compile rmats:

    ```bash
    conda activate rmats
    cd rmats-turbo/
    ./build_rmats
    ```

4. Modify GTF file path in the `SM_RNAseq.smk` file.
***Modify this with the actual name of the GTF file.

## 2. Prepare RNA-seq files

To correctly run the snakemake file, prepare each of the filenames inside the RNA-seq data in the follow way:

```text
GlobalSampleName_GroupName*_SampleNumber_1.fastq.gz
GlobalSampleName_GroupName*_SampleNumber_2.fastq.gz
```

*Group name should be annotated as `wt` for control samples, and is recommended to annotate the condition sample as `mut`.

### Example

```text
HEK293_wt_1_1.fastq.gz      HEK293_wt_1_2.fastq.gz
HEK293_mut_1_1.fastq.gz     HEK293_mut_1_2.fastq.gz
...
```

## 3. Modify Snakemake File paths

Edit `rMATS-pipeline.smk` variables:

- **raw_path**: For directory path to RNA-seq raw data to be used.
- **STAR_INDEX**: For directory path to STAR index.
- **GTF**: For directory path to genomic GTF reference data.

## 4. Run snakemake

Activate snakemake environments and run snakemake

```bash
conda activate snakemake
snakemake -c20 --use-conda -s rMATS-pipeline.smk -j4
```

## Preparing BAM files

- Done automatically by the snakemake
- Place BAM files in `bam_files` directory.
- Modify the two text files (`b1.txt` and `b2.txt`) with comma-separated path values to the BAM files from two different sample groups. If you are analyzing a single group, use only `b1.txt`.
- The format of the files should be modified as follows:

```bash
path/to/file_1.bam,path/to/file_2.bam,...,path/to/file_n.bam
```

**NOTE:** Very first run of the SnakeMake file ```SM_RNAseq.smk``` will initialize new conda environments. This may take a several minutes.
