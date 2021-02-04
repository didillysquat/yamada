# Yamada et al. 20201

## Overview
A pipeline for *de novo* transcriptome assembly using Trinity.

[Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used to get a QC overview of reads before processing and after processing.

Reads are trimmed with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

Error corrected with [Rcorrector](https://github.com/mourisl/Rcorrector).

Screened and filtered for rRNA using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [LSU+SSU](https://www.arb-silva.de/) SILVA databases.

Screened for taxonomy using [MMSeqs2](https://github.com/soedinglab/MMseqs2) in against the NCBI's taxonomically annotated nt database.

Assembled with [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki).

Assebly evautation statistics are produced using a collection of the perl scripts that come as part of the Trinity package and as documented [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment). As part of this evaluation, reads are mapped back to the assemblies using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and transcipts are quantified on a per sample basis to give raw and normalised read counts using [RSEM](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323) and the Trinity scripts. Metrics for number of genes, number of transcripts, percent of reads mapping back to the assemblies as well as ExN50 values are calculated amongst others.

[busco](https://busco.ezlab.org/) is also used to assess completeness of the assemblies using the eukaryota_odb10, alveolata_odb10 and stramenopiles_odb10 databases.

## Pipeline requirements
### Nextflow
To run this nextflow pipeline, [Nextflow](https://www.nextflow.io/) must be installed on your system.

### Docker
This pipeline uses Docker images/containers. Docker must therefore be installed on your system to run this pipeline. Docker images will automatically be downloaded if they are not found locally.

### NCBI's nt database with taxonomic annotations setup for MMSeqs
This pipeline runs taxonomic annotations using MMSeqs2 and the NCBI's nt database with taxonomic annotations. Specific setup guidelines are in the section "Create a seqTaxDB from an existing BLAST database" in the MMSeqs2 manual found [here](https://mmseqs.com/latest/userguide.pdf).

### SILVA LSU and SSU rRNA database bowtie2 indexed
This pipeline screens for rRNA reads and removes them. To do this you will need an MMSeq2 indexed SILVA database that contains the rRNA reads to be screened against. To build this database, concatenate the SILVA_138.1_LSUParc_tax_silva.fasta and SILVA_138.1_SSUParc_tax_silva.fasta files that can be downloaded from [here](https://www.arb-silva.de/no_cache/download/archive/release_138_1/Exports/). Then create the MMSeqs2 database and index.

## Running the pipeline
The pipeline is launched in the usual way and you will need to either modify the input paramaters in the `nextflow.config` file, or provide them on the command line in the usual way.
```
nextflow run ydenovo.nf
```
Or e.g.
```
nextflow run ydenovo.nf --raw_reads_dir </PATH/TO/READ/DIR/>
```
The pipeline is currently specific to the squenceing files of the Yamada et al. 2021 project. However, this can be modified to make the pipeline function on a generic set of RNA-seq input files.

N.B. due to a [bug](https://github.com/soedinglab/MMseqs2/issues/399) in MMSeq2, the evalue for taxonomic screening should not be more sensitive than '1e-35' (the current default)

## Outputs

The following output files are produced:

### fastqc_pre_trim

Contains the fastqc files of the input sequencing files before any preprocessing.

### fastqc_post_trim

Contains the fastqc files for after trimming with Trimmomatic.

### fastqc_post_correct

Contains the fastqc files for after error correction and removal of rRNA.

### cleaned_reads

Contains the fastq.gz read files after preprocessing.

### mmseqs_taxonomy_1e-XX

Contains the MMSeqs taxonomy outputfiles. Of particular intrest will be the .html files that give an interactive overview of the taxonomy results.

### trinity_assembly

Contains the Trinity assemblies and gene-to-transcript mapping files

### trinity_basic_stats

Contains files related to the basic statistics of the Trinity assemblies such as N50 and number of genes/transcripts.

### busco_stats

Contains the busco results.

### read_mapping_stats

Contains the results of mapping the reads back to the assemblies

### expression_quantification_rsem

Contains gene and transcript count data as well as ExN50 results.
