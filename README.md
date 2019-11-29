
Pamir: Discovery and Genotyping of Novel Sequence Insertions in Many Sequenced Individuals
======
Pamir detects and genotypes novel sequence insertions in single or multiple datasets of paired-end WGS (Whole Genome Sequencing) Illumina reads by jointly analyzing one-end anchored (OEA) and orphan reads.

# Table of contents
1. [Installation](#installation)
2. [Running Pamir](#running-pamir)
2. [Example](#example)
3. [Visualization](#visualization)
4. [Publications](#publications)
5. [Contact & Support](#contact-and-support)

# Installation

## Installation from Source
> *Prerequisite.* You will need g++ 5.2 and higher to compile the source code.

The first step to install Pamir is to download the source code from our 
[GitHub](https://github.com/vpc-ccg/pamir) repository. After downloading, 
change the current directory to the source directory ``pamir`` and run ``make`` and ``make install`` in
terminal to create the necessary binary files. 

```shell
git clone https://github.com/vpc-ccg/pamir.git --recursive
cd pamir
make
make install
``` 

# Running Pamir
## Prerequisites
Pamir's pipeline requires a number of external programs. You can either manually install them or 
take advantage of pamir's [conda](https://docs.conda.io/en/latest/) ``environment.yaml`` to 
install all the dependencies **except the assembler**:
```shell
conda env  create -f environment.yaml
source activate pamir-deps 
```

Dependencies | Version 
-------- |-----|
Python   | 3.x |
[samtools](http://www.htslib.org/) | >= 1.9 |
[mrsfast](https://github.com/sfu-compbio/mrsfast) | >= 3.4.0 |
[BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast\+/LATEST/) | >= 2.9.0+ |
[bedtools](https://bedtools.readthedocs.io/en/latest/) | >= 2.26.0 |
[bwa](https://github.com/lh3/bwa) | >= 0.7.17 |
[snakemake](https://snakemake.readthedocs.io/en/stable/) | >= 5.3.0 |
[RepeatMasker](http://www.repeatmasker.org/) | >= 4.0.9 |
[minia](https://github.com/GATB/minia) | >= 3.2.0 * |
[abyss](https://github.com/bcgsc/abyss) | >= 2.2.3 * |
[spades](https://github.com/ablab/spades) | >= 3.13 * |

***Note: You only need to install one of the assemblers.**

## Project Configuration
n order to run pamir, you need to create a project configuration file namely ``config.yaml``. 
This configuration consists of a number mandatory settings and some optional advance settings. 
Below is the list of the all the settings that you can set in your project.

|config-paramater-name | Type | Description|
|------------------------------|-----------|--------------------------------------------------------------------------------------------------------------------------------------|
| path                         | Mandatory | Full path to project directory.  |
| raw-data                     | Mandatory | Location of the input files (crams or bams) relative to ``path``.                                                         |
| population                   | Mandatory | Populuation/cohort name. Note that name cannot contain any space characters.                                                                                                                 |
| reference                    | Mandatory | Full path to the reference genome.                                                                                                        |
| read_length                  | Mandatory | Read length of the input reads.                                                                                                       |
| input                        | Mandatory | A list of input files per individual. Pamir 2.0 accepts BAM and CRAM files as input.       |
| analysis-base                | Optional  | Location of intermediate files relative to ``path``. default: ``{path}/analysis``|
| results-base                 | Optional  | Location of final results relative to the ``path``. default: ``{path}/results``  |
| min_contig_len               | Optional  | Minimum contig length from the external assembler to use. default: ``read_length``  |
| assembler                    | Optional  | External assembler to use (``minia``, ``abyss``, ``spades``) default: ``minia``                                                                      |
| assembler_k                  | Optional  | kmer to use for external assembler. default: 47                                                                                      |
| pamir_partitition_per_thread | Optional  | Number of internal pamir jobs to be completed per thread. This is an advanced settings, modifying this can heavily affect the performance. Too small or too large may affect the performance negatively.  default: 1000                     |
| blastdb                      | Optional  | Full path to blast database to remove possible contaminants from the data.  |
| centromeres | Optional | Full path to the file in bed format that contains centromeres locations. The calls in these regions will not be reported |
| align_threads                | Optional  | number of threads to use for alignment jobs. default: 16                                                                             |
| assembly_threads             | Optional  | number of threads to use for assembly jobs. default: 62                                                                              |
| other_threads                | Optional  | number of threads to use for other jobs. default: 16                                                                                 |
| minia_min_abundance          | Optional  | minia's internal assembly parameter. default: 5|



The following a an example of ``config-yaml`` with two individuals. 
```yaml
path:
    /full/path/to/project-directory
raw-data:
    raw-data
read_length:
    100
reference:
    /full/path/to/the/reference.fa
population:
    my-pop
input:
 "samplename1":
  - A.cram
 "samplename2":
  - B.cram/bam
```

Now, to run pamir on such a config file, you have to run the following command.
```shell
pamir.sh  --configfile /path/to/config.yaml
```

Since, ``pamir.sh`` is internally utilizing ``snakemake``, you can pass any additionak ``snakemake`` parameters to ``pamir.sh``. Here are some examples:
```shell
pamir.sh  --configfile /path/to/config.yaml -j [number of threads] 
pamir.sh  --configfile /path/to/config.yaml -np [Dry Run] 
pamir.sh  --configfile /path/to/config.yaml --forceall [rerun all steps regardless of the current stage]
```

### Read Length
Now Pamir only accepts WGS datasets in which two mates of all reads are of **equal length**.

## Output Formats
Pamir will generate the following structure.
Pamir generates a [VCF file](https://samtools.github.io/hts-specs/VCFv4.2.pdf) for detected novel sequence insertions. 

```
[path]/
├── raw-data                       -> OR [raw-data]
│   ├── A.cram
│   ├── B.cram
├── analysis                       -> OR [analysis-base]
│   └── my-pop
└── results                        -> OR [results-base]
    └── my-pop
        ├── index.html             -> Summary fo events
        ├── summary.js             -> Summary required by index.html
        ├── data.js                -> Data required by index.html
        ├── events.repeat.bed      -> annotation of repeats for detected eveents
        ├── events.fa              -> all the detected events with 1000bp flanking region
        ├── events.fa.fai          -> index of events.fa
        └── ind
            ├── A
            │   ├── events.bam     -> mapping of the reads in the events region
            │   ├── events.bam.bai -> index
            │   ├── events.bed     -> location of events
            │   └── events.vcf     -> genotyped insertion calls
            ├── B
            │   ├── events.bam
            │   ├── events.bam.bai
            │   ├── events.bed
            │   └── events.vcf
```

# Example
```shell
curl -L https://ndownloader.figshare.com/files/20021672?private_link=ce133e9aaa4d38a6d6b0 --output example.tar.gz
tar xzvf example.tar.gz
cd example; ./configure.sh
pamir.sh -j16 --configfile config.yaml
```

# Visualization
``index.html`` provides a quick way of looking at general overview of events. It is an alternative to working with vcf files in a friendly fashion.
If you start your IGV, you can easily jump back and forth investigating your events from ``index.html``.  


# Publications
**Discovery and genotyping of novel sequence insertions in many sequenced individuals.** P. Kavak*, Y-Y. Lin*, I. Numanagić, H. Asghari, T. Güngör, C. Alkan‡, F. Hach‡. [***Bioinformatics*** (ISMB-ECCB 2017 issue), 33 (14): i161-i169, 2017.](https://doi.org/10.1093/bioinformatics/btx254)

# Contact and Support

Feel free to drop any inquiry at the [issue page](https://github.com/vpc-ccg/pamir/issues)    .

