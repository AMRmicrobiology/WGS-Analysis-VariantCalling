# WGS-Analysis-VariantCalling
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![license-shield]][license-url]

## Introduction
This repository hosts an advanced pipeline build with Nextflow for whole-genome sequencing (WGS) analysis and genetic variant calling, specifically optimized for Illumina sequencing data of bacterial genomes. It is designed to offer an automated, reproducible, and scalable solution for processing large-scale genomic data in clinical microbiology research.

## Contents
- [Installation](#installation)
- [Process](#process)
    - [de-novo](#de-novo)
    - [refence-genome](#reference-genome)
- [Requirements](#requirements)
- [How to Use It](#how-to-use-it)
    - [Parameters](#parameters)
- [Project Structure]()
- [References](#reference)

## Installation
Set Up the Environment:
- Install Nextflow: [Installation Guide](https://github.com/nextflow-io/nextflow)
- Install [Docker](https://github.com/docker/docker-install) or [Singularity](https://github.com/sylabs/singularity-admindocs/blob/main/installation.rst) for container support.
- Ensure Java 8 or higher is installed.

## Process:
The pipeline includes the following steps:

### Quality Control (Q.C.)

1. Quality Control: Assessment of raw sequencing data using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to evaluate read quality. Removal of low-quality bases and adapter sequences with [Trimmomatic](https://github.com/usadellab/Trimmomatic).

### de-novo

2. de novo assembly using [spades]().
3. Assembly evaluation with [QUAST](https://bioinf.spbau.ru/quast), structural quality metrics of the assembly and [BUSCO](https://github.com/metashot/busco), evaluation of biological completeness.
4. Alignment using [BWA](https://github.com/bwa-mem2/bwa-mem2).


### reference genome

2. Alignment using [BWA-MEM](https://github.com/bwa-mem2/bwa-mem2), [samtools](https://github.com/samtools/samtools) 
3. Q.Control
4. [MultiQC](https://github.com/MultiQC/MultiQC)


### Process to indentify and filter variants:
  5.1.  Anotation using [Prokka](https://github.com/tseemann/prokka), [Bakta](https://github.com/oschwengers/bakta)

  5.2 . Variant Identification: Detection of single nucleotide polymorphisms (SNPs) and insertions/deletions (indels) using [PicardTools](), [GATK](https://github.com/broadinstitute/gatk) or [FreeBayes](https://github.com/freebayes/freebayes).

  5.3. Variant Filtering: Application of quality filters to obtain high-confidence variant calls.
6. Post-Alignment Analysis:

  6.1. Genetic variant annotation usign [SnpEff](http://pcingola.github.io/SnpEff/), a toolbox for annotating and predicting the functional effects of genetic variants on genes and proteins.

  6.2. Mass screening of contigs for antimicrobial resistance or virulence genes using [ABRIcate](https://github.com/tseemann/abricate).

  6.3. Identification of acquired antimicrobial resistance genes and point mutations in protein and/or assembled nucleotide sequences using [AMRFinder](https://github.com/ncbi/amr).

  6.4. Predict of Antibiotic Resistance Genes using [DeepARG](https://github.com/gaarangoa/deeparg).

  6.5. 

## Requirements


## How to use it?

Clone the Repository:

```
git clone https://github.com/AMRmicrobiology/WGS-Analysis-VariantCalling.git
cd WGS-Analysis-VariantCalling
```
Run the pipeline using the following command, adjusting the parameters as needed:

CLINICAL

```
nextflow run main.nf --mode clinical --input '/path/to/data/*.fastq.gz' --outdir './out' -profile docker
```

REFERENCE STRAIN

```
nextflow run main.nf --mode refrence --input '/path/to/data/*.fastq.gz' --personal_ref '/path/to/bacterial_genome.fasta' --outdir './out' -profile docker
```

### Parameters

--mode: dpende of the analysis clinical/reference.

--input: Path to input FASTQ files generated by Illumina sequencing (file format: .fastq.gz).

--personal_ref: Path to the bacterial reference genome FASTA file.

--outdir: Directory where the results will be stored (default: out).

-profile: Specifies the execution profile (docker, singularity, standard, HPC).

add:

--workDir: Path to the temporary work directory where files will be stored (default: ./work).



## Project Structure


[contributors-shield]: https://img.shields.io/github/contributors/jimmlucas/DIvergenceTimes.svg?style=for-the-badge
[contributors-url]: https://github.com/jimmlucas/DIvergenceTimes/graphs/contributors

[forks-shield]: https://img.shields.io/github/forks/jimmlucas/DIvergenceTimes.svg?style=for-the-badge
[forks-url]: https://github.com/jimmlucas/DIvergenceTimes/network/members

[stars-shield]: https://img.shields.io/github/stars/jimmlucas/DIvergenceTimes.svg?style=for-the-badge
[stars-url]: https://github.com/gjimmlucas/DIvergenceTimes/stargazers

[issues-shield]: https://img.shields.io/github/issues/jimmlucas/DIvergenceTimes.svg?style=for-the-badge
[issues-url]: https://github.com/jimmlucas/DIvergenceTimes/issues

[license-shield]: https://img.shields.io/github/license/jimmlucas/DIvergenceTimes.svg?style=for-the-badge
[license-url]: https://github.com/jimmlucas/DIvergenceTimes/blob/master/LICENSE.txt

## Reference:

[In Silico Evaluation of Variant Calling Methods for Bacterial Whole-Genome Sequencing Assays](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10446864/)

[Recommendations for clinical interpretation of variants found in non-coding regions of the genome](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9295495/)












