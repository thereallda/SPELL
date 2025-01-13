# Description

The `ChIPseq_Pipeline.sh` script integrates the common steps of ChIP-seq analysis such as 

1. QC(adapters trimming included)

2. Mapping

3. sam2bam, sort and index

4. bam index

 while the `run_bowtie2_pipeline.sh` can handle multi-samples pipeline analysis.

 

# Dependences

```shell
fastqc #v0.11.8
cutadapt #v2.10
bowtie2 #2.3.5
samtools #1.9
```



## Usage

```shell
  bash BOWTIE2_Pipeline.sh [options] --fastq1 <fastq1_file>  -d <project_dir> -o <output_dir> -a <forward_strand 3 adapter> --ref <reference_genome>
```



 ### Options

  `--fastq1`	fastq1 file
  `-a`	Forward_adapt sequence
  `-d`	Project directory
  `-o`	Output directory to store the results
 ` --ref`	Index using for bowtie2 mapping

Options:
`  -t, --thread`	Numbers of threads, defult: 1
`--min-length`	 Discard reads shorter than min_length, defult:30
  `--pair`	 Flag to turn on pair-end mode
  `--multi`	Flag to turn on multi-fastq/samples mode
  `-h`	 show help message

Pair-end options:
  `--fastq2`	fastq2 file when --pair is on
  `-A`	Reverse_adapt sequence



## Multi-samples pipeline usage

```shell
bash run_bowtie2_pipeline.sh -d <project_dir> -o <output_dir> -a <forward_strand 3 adapter> --ref <reference_genome> --samples <Dir_names_of_samples>
```

- Output
  1. The trimmed fastqs in `<project_dir>/data/clean` following the structure of your fastq directory

  ```shell
  data/clean/
  ...
  ```
  
  2. The QC results of raw reads and trimmed reads in ` <output_dir>/QC/{clean,raw}`
  
  ```shell
  results/QC/
  |-- clean
  ...
  `-- raw
  ...
  ```

  3.  The bowtie2 mapping and samtools results in `<output_dir>/bowtie/`
  
  ```shell
  results/bowtie/
  
  ...
  
  ```



## Test


I tested the scripts by using  the ChIP-seq data in project `Histone_modification_fly_development`

```shell
## Directory strcture
lidean@admin 21:38 ~/LNlab_project/test
$tree data/fastq/
data/fastq/
|-- H3K27me3-Repeat2
|   |-- 20.fastq.gz
|   `-- 2.fastq.gz
|-- H3K4me3-Repeat1
|   |-- 3-20.fastq.gz
|   `-- 3-2.fastq.gz
`-- H3K4me3-Repeat4
    |-- N20.fastq.gz
    `-- N2.fastq.gz
```



```shell
bash run_bowtie2_pipeline.sh -d ~/LNlab_project/test -o ~/LNlab_project/test/results -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --ref /public/genomes/bowtie2index/dm6/genome -t 8 --samples H3K27me3-Repeat2,H3K4me3-Repeat1,H3K4me3-Repeat4
```



> DEVELOP LOG
>
> --version 1.0
>
> 2020-08-19
>
> By now, the single-end mode of `BOWTIE2_Pipeline.sh` and its multi-samples wrapper  `run_bowtie2_pipeline.sh` have been run successfully, while the pair-end mode have not been tested. 
>
> Using `--multi` flag to turn on multi-samples mode in `BOWTIE2_Pipeline` and depending on  `${sample}`  variable exported by `run_bowtie2_pipeline.sh` with no option in  `BOWTIE2_Pipeline`  handling it which are somehow inconsistency and needed to be fixed.
>
> In addition, the for loop now just call the script in a serial mode rather than in a parallel mode.
>
> --version 1.1
>
> 2020-08-25
>
> In the development of RNA-seq pipeline scripts, I fixed the pair-end mode problems by adding `if condition` in the `run_pipeline.sh` .



## v1.2 

### Update

- In this version, the pipeline scripts only perform QC on trimmed data for the reason that the QC of raw data seems to be redundant for automatic pipeline. 
- Also, the default setting of adapter options are both Illumina universal adapter  `AGATCGGAAGAG`
- The `multi` option is now removed
- Add `-i` parameter to specify the data directory that holds raw fastqs and trimmed fastqs

For multiple fastqs pipeline wrapper `run_star_pipeline.sh`, the `samples` option is also removed. The scripts use for loop to pass `fastqs` and only manage the results by fastq's ID without sample information.



### Example

```shell
bash run_bowtie2_pipeline.sh -d ~/LNlab_project/ChIP-seq_jiang_20210408 -i ~/LNlab_project/ChIP-seq_jiang_20210408/data -o ~/LNlab_project/ChIP-seq_jiang_20210408/results --ref /public/Reference/human/index/bowtie2_index/GRCh38.primary_assembly/GRCh38 -t 48
```






## Authors

- Dean Li

## Version

- 1.0
- 1.1
  - pair-end mode tested with PE RNA-seq

- 1.2
  - Remove "Step01 fastqc to raw data" in `STAR_pipeline.sh`
  - Remove `--multi` option  in `STAR_pipeline.sh`
  - Change `-a, -A` adapter options to default: `AGATCGGAAGAG`
  - Remove `--sample` option  in `run_star_pipeline.sh`
  - `-i` 