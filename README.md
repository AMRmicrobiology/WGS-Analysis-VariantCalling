# WGS-Analysis-VariantCalling
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![license-shield]][license-url]

This repository hosts an advanced pipeline developed with Nextflow for whole-genome sequencing (WGS) analysis and genetic variant calling, specifically optimized for Illumina sequencing data of bacterial genomes. It is designed to provide an automated, reproducible, and scalable solution for processing large-scale genomic data in microbiology research.
## Table of Contents
- [Process](#process)
- [Requirements](#requirements)
- [How to Use It](#how-to-use-it)
-- [Parameters](#parameters)
- [Project-Structure]()
- [References](#reference)

## Process:
The pipeline includes the following steps:

1. Quality Control: Assessment of raw sequencing data using FastQC to evaluate read quality. Removal of low-quality bases and adapter sequences with Trimmomatic. 
2. Alignment using BWA-MEM.
4. Post-Alignment Processing
5. Variant Calling: Detection of single nucleotide polymorphisms (SNPs) and insertions/deletions (indels) using GATK or FreeBayes.
6. Variant Filtering: Application of quality filters to obtain high-confidence variant calls.
7. Annotation: Functional annotation of variants with SnpEff to predict their impact on genes and proteins.

## Requirements
Set Up the Environment:
- Install Nextflow: [Installation Guide](https://github.com/nextflow-io/nextflow)
- Install [Docker](https://github.com/docker/docker-install) or [Singularity](https://github.com/sylabs/singularity-admindocs/blob/main/installation.rst) for container support.
- Ensure Java 8 or higher is installed.

## How to use it?

Clone the Repository:

```
git clone https://github.com/AMRmicrobiology/WGS-Analysis-VariantCalling.git
cd WGS-Analysis-VariantCalling
```
Run the pipeline using the following command, adjusting the parameters as needed:

```
nextflow run main_clinical.nf \
  --reads '/path/to/data/*.fastq.gz' \
  --genome '/path/to/bacterial_genome.fasta' \
  --outdir './out' \
  -profile docker
```

### Parameters

main: [main_clinal.nf](main_clinial.nf) includes the pipeline for de novo assembly; [main_reference.nf](main_reference.nf) includes the pipeline for analysis using a reference genome.

--reads: Path to input FASTQ files generated by Illumina sequencing (file format: .fastq.gz).

--genome: Path to the bacterial reference genome FASTA file.

--outdir: Directory where the results will be stored (default: out).

-profile: Specifies the execution profile (docker, singularity, standard, HPC).

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












