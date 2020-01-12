# scion
### Single Cell Bioinformatics

This pipeline is a wrapper of various tools associated with single cell RNA-seq analysis. The pipeline is optimized for 3' based library preparation methods in NGS for single cell. This pipeline is specifically optimized for CEL-Seq2 libraries that have pseudo-paired reads.

Current library prep requirements are:
  * CEL-Seq2
  * Already lane merged fastq files that result in a single pair of fastqs.
    * Read 1 FASTQ
      * The UMI in position 1 - 6
      * The cell barcode in position 7-12
    * Read 2 FASTQ
      * The transcript 1-50 (or the max length of the second record)
  * Sequenced on Illumina platform
  * Relatively low throughput due to sample sheet requirements
      
## Getting Started

The pipeline is dockerized and contains all the necessary software dependencies to run quality control on raw reads, perform alignment and UMI counting, and run a basic differential expression analysis with a simple experiment design.

### Prerequisites

* Refer to the Sample Sheet README Wiki for creating the Sample Sheet

* Docker
  * https://docs.docker.com/v17.12/install/
  * Install per their instructions
  * If installing on Linux machine, must have superuser access to create and add users as a group
    * https://docs.docker.com/install/linux/linux-postinstall/
  * Docker commands in this repository are **not** ran as sudo
  * If running this repository on AWS, docker may fail due to AWS intricacies
  
### Host machine directory setup

On the local machine:
  * Create a directory that contains the following:
      * Reference_Genome_FASTA file
      * Reference_Annotation gtf file
      * Compressed and Lane Merged Read 1 FASTQ (gz)
      * Compressed and Lane Merged Read 2 FASTQ (gz)
      * SampleSheet.csv
  * Refer to this repository's wiki **Sample Sheet README** on creating the sample sheet
  * Sample sheet **must** be named SampleSheet.csv

After creating directory and adding the files, change into the directory and clone this repository.
```
cd pipeline_example
git clone https://github.com/tbeaux18/sc
```
This results in the following directory structure.
```
root:pipeline_example example$
.
├── scion
├── SampleInfo.csv
├── ExperimentDesign.csv
├── read1.fastq.gz
├── read2.fastq.gz
├── ref_genome.fa
└── ref_genome_annotation.gtf
```

To run the pipeline:
```
cd scion
./run_pipeline.sh
```

### Pipeline Procedure
By running this command:
  * Docker images will build
  * Docker runs the **ubuntur35:pipeline** container
  * Docker mounts the parent directory of scion
  * Sample Sheet is parsed and exchanges information with necessary input files such as zUMIs.yaml file
  * Quality Control of Raw Data
  * STAR Index is built with the provided FASTA, the GTF is not included in the index build
  * zUMIs pipeline begins
    * Stages
      * Filtering
      * Mapping
      * Counting
      * Summarising
    * Output directory
      * <sample_sheet_basename>_zumi_output
        * BAM Files
        * zUMI log files
        * Count matrices
        * Alignment plots


### Limitations and Future Development

* Requires zUMIs V2.4.0
  * Currently UMI collapsing is done using 1 core; slows down counting immensley.
  * Large amounts of cells required to multi-thread
* Robust error handling will be included
* Sample sheet and pipeline usage will expand to include:
  * Different library preparation methods
  * More zUMIs parameters for customization
  * Different entrypoints


### Code languages to build repository

* [python3] (https://www.python.org/download/releases/3.0/)
* [bash] (https://www.gnu.org/software/bash/)
* [R] (https://www.r-project.org/)
* [Docker] (https://www.docker.com/)

### Docker container software

* [FASTQC] (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Cutadapt] (https://cutadapt.readthedocs.io/en/stable/) - python3
* [zUMIs] (https://github.com/sdparekh/zUMIs) - V2.4.0
* [STAR] (https://github.com/alexdobin/STAR) - V2.7
* [Rsubread] (https://bioconductor.org/packages/release/bioc/html/Rsubread.html) - featureCounts

## Authors

* **Timothy Baker**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Dr. Jennifer Beshel, PhD
* Dr. Catherine Putonti, PhD
* Loyola University Chicago
* CEL-Seq2 authors


