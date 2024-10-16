# WGS-Analysis-VariantCalling
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![license-shield]][license-url]

## Introduction
This repository hosts an advanced pipeline build with Nextflow for whole-genome sequencing (WGS) analysis and genetic variant calling, specifically optimized for Illumina sequencing data of bacterial genomes. It is designed to offer an automated, reproducible, and scalable solution for processing large-scale genomic data in clinical microbiology research.

## Contents
- [Pipeline summary](#pipeline-summary)
    - [*de-novo*](#de-novo)
    - [Refence genome](#reference-genome)
- [Installation](#installation)
- [How to Use It](#how-to-use-it)
    - [Parameters](#parameters)
- [References](#reference)

## Pipeline summary:
The pipeline includes the following steps:

1. **Quality Control**: Assessment of raw sequencing data using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to evaluate read quality. Removal of low-quality bases and adapter sequences with [FastP](https://github.com/OpenGene/fastp) followed again by [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://github.com/MultiQC/MultiQC) to summarise the input data.

At this point, the two modes available in the pipeline differ on the input reference genome. You can perform the variant calling using a [*de novo*](#de-novo) assembled reference strains or an already available [reference genome](#reference-genome). 

-  #### *De-novo*

    - **Assembly**: After quality control as previously described, *de novo* assembly using [SPAdes](https://github.com/ablab/spades).
    -  **Quality assembly assessment**: Structural quality metrics of the assembly using [QUAST](https://bioinf.spbau.ru/quast) and evaluation of biological completeness with [BUSCO](https://github.com/metashot/busco).
    -   **Anotation**: Genome anotation using [Prokka](https://github.com/tseemann/prokka) and [Bakta](https://github.com/oschwengers/bakta).

-  #### Reference genome
    >The pipeline includes an script to download the reference genome.
    ```
    bash ./workflow/bin/download_reads.sh
    ```

After inputed the reference genome, the pipeline follows the same steps for both modes:

2. **Alignment**: Alignment against the selected reference genome with [BWA-MEM](https://github.com/bwa-mem2/bwa-mem2) and [samtools](https://github.com/samtools/samtools).
3. **Quality control**: Alignment quality control using [QUAST](https://bioinf.spbau.ru/quast).
4. **Aggregation of quality reports**: [MultiQC](https://github.com/MultiQC/MultiQC)

5. **Variant calling and filtering**:

    -  **Variant Identification**: Detection of single nucleotide polymorphisms (SNPs) and insertions/deletions (indels) using [PicardTools](https://broadinstitute.github.io/picard/), [GATK](https://github.com/broadinstitute/gatk) and/or [FreeBayes](https://github.com/freebayes/freebayes).

    -  **Variant Filtering**: Application of quality filters to obtain high-confidence variant calls ([*see Parameters*](#parameters)).

    -  **Genetic variant annotation**: Using [SnpEff](http://pcingola.github.io/SnpEff/), a toolbox for annotating and predicting the functional effects of genetic variants on genes and proteins.

7. **Post-Alignment Analysis**:
    
    - Mass screening of contigs for antimicrobial resistance or virulence genes using [ABRIcate](https://github.com/tseemann/abricate).

    -  Identification of antimicrobial resistance genes and point mutations in protein and/or assembled nucleotide sequences using [AMRFinder](https://github.com/ncbi/amr).

    - Prediction of Antibiotic Resistance Genes using [DeepARG](https://github.com/gaarangoa/deeparg).
 

## Installation
The prerequisites to run the pipeline are:
- Install [Nextflow](https://github.com/nextflow-io/nextflow)
- Install [Docker](https://github.com/docker/docker-install) or [Singularity](https://github.com/sylabs/singularity-admindocs/blob/main/installation.rst) for container support
- Ensure [Java 8](https://github.com/winterbe/java8-tutorial) or higher is installed

Clone the Repository:

```
# Clone the workflow repository
git clone https://github.com/AMRmicrobiology/WGS-Analysis-VariantCalling.git

# Move in it
cd WGS-Analysis-VariantCalling
```
<!-- compl -->
### Local (conda)

  ```
  conda create -n bacteriano -f enviromentWGS.yaml
  conda activate bacteriano
  ```



## How to use it?

Run the pipeline using the following command, adjusting the parameters as needed:

*DE NOVO*
>[!IMPORTANT]
The name of the paired-end reads of the reference sample must be labelled as **1** (e.g. AB**1**_1 / AB**1**_2)

```
nextflow run main.nf --mode novo --input '/path/to/data/*.fastq.gz' --genome_name_db ¨Acinetobacter_baumanii_clinical¨ --outdir './out' -profile <docker/singularity/conda>
```

*REFERENCE GENOME*
```
nextflow run main.nf --mode refrence --input '/path/to/data/*.fastq.gz' --personal_ref '/path/to/bacterial_genome.fasta' --outdir './out' -profile <docker/singularity/conda>
```

### Parameters

--mode: Depends on the analysis novo/reference/conda.

--input: Path to input FASTQ paired-end files generated by Illumina sequencing (file format: .fastq.gz).

--outdir: Directory where the results will be stored (default: out).

-profile: Specifies the execution profile (docker, singularity or local).

--genome_name_db (only for --mode novo): Name of the organism that will name the databse in SnpEFF.

--personal_ref (only for --mode reference): Path to the bacterial reference genome FASTA file.


#### Optional parameters

--workDir: Path to the temporary work directory where files will be stored (default: ./work).

##### Trimming

--cut_front (15):

--cut_tail (20):

--cut_mean_quality (20):

--length_required (50):

##### Filter

--qual_snp ("QUAL < 50.0 || MQ < 25.0 || DP < 30"): gattaka

--qual_indel ("QUAL < 200.0 || MQ < 25.0 || DP < 30"):


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

[An ANI gap within bacterial species that advances the definitions of intra-species units](https://journals.asm.org/doi/10.1128/mbio.02696-23)
<!-- ADD REFERENCES -->











