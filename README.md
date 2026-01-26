# AltSS

Alternative Splicing Snakemake pipeline, to process data from RNA-Seq to alternative splicing prediction. Only working with rMATS [1] and SplAdder [2] tools so far.

[1]<https://github.com/Xinglab/rmats-turbo><br>
[2]<https://spladder.readthedocs.io/en/latest/index.html>

## Requirements

Tested in WSL2 and linux server

- Conda 24.11.3
- GTF and STAR index files for respective genomic release (Genecode or Ensembl).

## 1.1. Prepare rMATS environment

1. Download `rmats-turbo` in a folder with the same name.

    ```bash
    git clone https://github.com/Xinglab/rmats-turbo.git
    ```

2. Create conda environment for Snakemake and rMATS.

    ```bash
    conda env create -n snakemake -f envs/snakemake.yml
    conda env create -n rmats -f envs/rmats.yml
    ```

3. Activate rmats environment, and compile rmats:

    ```bash
    conda activate rmats
    cd rmats-turbo/
    ./build_rmats
    ```

## 1.2. Prepare SplAdder environment

1. Create conda environment for Snakemake and SplAdder.

    ```bash
    conda env create -n snakemake -f envs/snakemake.yml
    conda env create -n spladder -f envs/spladder.yml
    ```

## 2. Prepare RNA-seq files

To correctly run the snakemake file, prepare each of the filenames inside the RNA-seq data in the follow way:

### For paired reads

```text
GlobalSampleName_GroupName*_SampleNumber_1.fastq.gz
GlobalSampleName_GroupName*_SampleNumber_2.fastq.gz
```

### For single reads

```text
GlobalSampleName_GroupName*_SampleNumber.fastq.gz
```

*Group name corresponds to the condition and control samples.

### Example

```text
HEK293_wt_1_1.fastq.gz      HEK293_wt_1_2.fastq.gz
HEK293_mut_1_1.fastq.gz     HEK293_mut_1_2.fastq.gz
...
```

## 3. Configure Pipeline

Edit `config.yaml` with the following variables:

### Required Parameters

- **raw_path**: Directory path to RNA-seq raw data to be used.
- **gtf**: Path to genomic GTF reference file.
- **fasta**: Path to genome FASTA file (required for generating a new STAR index).
- **outdir**: Output directory for results.
- **read_length**: Read length of your sequencing data. Default: 50.

### Optional Parameters

- **read_type**: Either "paired" or "single". Default: "paired".
- **control_name**: Name of the control group for identification of the group in the name file
- **nthread**: Number of threads to use. Default: 8 cores.

### STAR Index Generation

The pipeline can automatically generate a STAR index. This requires significant RAM (typically 30-40GB for human genome). Ensure your system has sufficient memory or use lower number of cores (nthread).

## 4. Run snakemake

Activate snakemake environments and run snakemake

```bash
conda activate snakemake
snakemake -c --use-conda -s rMATS-pipeline.smk --configfile config.yaml
```

Or for SplAdder:

```bash
conda activate snakemake
snakemake -c --use-conda -s SplAdder-pipeline.smk --configfile config.yaml
```

**NOTE:** Very first run of the SnakeMake file, it will initialize new conda environments. This may take a several minutes.

## Troubleshooting

### ERROR: libgsl.so.25: cannot open shared object file: No such file or directory

If rMATS finishes with the aforementioned error, it means that the program is not locating the correct path for that function. Modify this by using the follow bash command, while in the snakemake environment.

```bash
export LD_LIBRARY_PATH=path/to/rmats_pipeline/.snakemake/conda/environement_created_for_rmats/lib
```
