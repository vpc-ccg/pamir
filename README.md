
Pamir: Discovery and Genotyping of Novel Sequence Insertions in Many Sequenced Individuals
======
Pamir detects and genotypes novel sequence insertions in single or multiple datasets of paired-end WGS (Whole Genome Sequencing) Illumina reads by jointly analyzing one-end anchored (OEA) and orphan reads.

# Table of contents
1. [Installation](#installation)
2. [Output Formats](#output-formats)
3. [Project Configuration](#project-configuration)
4. [Publications](#publications)
5. [Contact & Support](#contact-and-support)

## Installation
Source code of Pamir can be downloaded from [GitHub](https://github.com/vpc-ccg/pamir). To begin with, you will need to set up several external tools as described below.

### Prerequisite
Pamir relies on specific version of the following tools:
1. g++ >= 5.2.0
2. Python3
3. [samtools](http://www.htslib.org/) >= 1.9
4. [mrsfast](https://github.com/sfu-compbio/mrsfast) >= 3.4.0
5. [BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast\+/LATEST/) >= 2.9.0+
6. [bedtools](https://bedtools.readthedocs.io/en/latest/) >= 2.26.0
7. [bwa](https://github.com/lh3/bwa) >= 0.7.17
8. [snakemake](https://snakemake.readthedocs.io/en/stable/) >= 5.3.0
9. Assembler
* [minia](https://github.com/GATB/minia) >= 3.2.0 Or
* [abyss](https://github.com/bcgsc/abyss) >= 2.2.3 Or
* [spades](https://github.com/ablab/spades) >= 3.13.1

#### Or you can use environment.yaml with conda (Except Assembler)
    conda env  create -f environment.yaml
    source activate pamir-deps 




   You also need to download the latest BLAST nt database to /your/path/to/ncbi-blast-2.5.0+/db/ (see *Compilation and Configuration* below) for contamination detection. 


```
mkdir /dir/to/blast/db
cd /dir/to/blast/db
/path/to/blast/bin/update_blastdb.pl nt
```

### Compilation and Configuration
To install Pamir, you need to first fetch Pamir from our [git repository](https://github.com/vpc-ccg/pamir) or download the corresponding compressed files. 
```
git clone https://github.com/vpc-ccg/pamir.git
cd pamir
```

Make the project and move running script and executable folder somewhere in your PATH.
```
pamir$ make -j
pamir$ mv pamir /usr/bin/
pamir$ mv pamir.sh /usr/bin/
```

## Project Configuration
Pamir is designed to detect novel sequence insertions based on one-end anchors (OEA) and orphans from paired-end Whole Genome Sequencing (WGS) reads. A .yaml configuration file with the following fields.

1. path: full path to your project working directory
2. raw-data: Where are the crams(bams) are stored relative to the path(1)
3. analysis-base: Where intermediate files should be stored relative to the path(1)
4. results-base: Where final results should be stored relative to the path(1)
5. population: population name/id
6. reference: Whole genome reference path
7. min\_contig\_len: Minimum contig length from the external assembler to use (Should be same as read length)
8. assembler: external assembler to use (minia, abyss, spades)
9. assembler\_k: kmer to use for external assembler
10. assembly\_segmenter\_count: Number of processes to split the pamir partition work (Higher = more parallelism + more overhead) 
11. blastdb: Path to blast database
12. read\_length: read length of the input
13. align\_threads: number of threads to use for alignment jobs
14. assembly\_threads: number of threads to use for assembly jobs
15. other\_threads: number of threads to use for other jobs
16. input: a list of input files per individual like following, where files should be in the (path(1))/(raw-data(4))/. This version of pamir accepts BAM and CRAM files as input.

```
input:
  "A":
   - A.cram
  "B":
   - B.cram
  "C":
   - C.cram
```

Pamir can be run with the following commands, this version of pamir requires 2 steps to be run in order. Where first command create the partition files, and second command genotypes the insertions.
```
pamir.sh  --configfile /path/to/config.yaml   [path]/[analysis-base]/[population]/pamir/partition/done
pamir.sh  --configfile /path/to/config.yaml
```

You can pass any snakemake parameter to pamir.sh.
```
-j [number of threads]
-np [Dry Run]
--forceall [rerun all steps regardless of the current stage]
etc.
```


#### Read Length
Now Pamir only accepts WGS datasets in which two mates of all reads are of **equal length**.


## Output Formats
Pamir generates a [VCF file](https://samtools.github.io/hts-specs/VCFv4.2.pdf) for detected novel sequence insertions. You can run genotyping for each sample after obtaining the VCF file by:

Your project will be in the following structure

[path] ->  | -> [raw-data] -> |-> A.cram
           |                  |-> B.cram
           |                  |-> C.cram
           |
           | -> [analysis-base] ->  [population] -> [...Intermediate Files...]
           |
           | -> [results-base]  -> [population] -> |-> ind  ------>|->A ----|->events.bam
                                                   |-> events.ref  |->B ..  |->events.vcf                         
                                                                   |->C ..  |->events.bed

events.vcf contains genotyped insertions calls.
events.fa  is the reference build from called insertions + flanking regions from the genome
events.bam is the reads from flanking regions of the calls
events.bed is the bed file showing the insertion positions on events.fa


## Small working example
    curl -L https://ndownloader.figshare.com/files/18706463?private_link=42900675d70a9a2282e8 --output small-example.tar.gz
    tar xzvf small-example.tar.gz
    cd small-example; ./configure.sh
    pamir.sh -j16 --configfile config.yaml $(pwd)/analysis/sim/pamir/partition/done
    pamir.sh -j16 --configfile config.yaml

## Publications
**Discovery and genotyping of novel sequence insertions in many sequenced individuals.** P. Kavak*, Y-Y. Lin*, I. Numanagić, H. Asghari, T. Güngör, C. Alkan‡, F. Hach‡. [***Bioinformatics*** (ISMB-ECCB 2017 issue), 33 (14): i161-i169, 2017.](https://doi.org/10.1093/bioinformatics/btx254)

## Contact and Support

Feel free to drop any inquiry at the [issue page](https://github.com/vpc-ccg/pamir/issues)    .

