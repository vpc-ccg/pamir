Pamir: Discovery and Genotyping of Novel Sequence Insertions in Many Sequenced Individuals
======
Pamir detects and genotypes novel sequence insertions in single or multiple datasets of paired-end WGS (Whole Genome Sequencing) Illumina reads by jointly analyzing one-end anchored (OEA) and orphan reads.


## Installation
Source code of Pamir can be downloaded from [GitHub](https://github.com/vpc-ccg/pamir). To begin with, you will need to set up several external tools as described below.

### Prerequisite
Pamir relies on specific version of the following tools:
1. g++ 4.9.0 or higher
2. Python 2.7 or higher (for the package *argparse*)
4. [Velvet](https://github.com/dzerbino/velvet) 1.2.10 or higher
5. [BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 2.3.0+ or higher

   You also need to download the latest BLAST nt database to /your/path/to/ncbi-blast-2.5.0+/db/ (see *Compilation and Configuration* below) for contamination detection. 

```
mkdir /dir/to/blast/db
cd /dir/to/blast/db
/path/to/blast/bin/update_blastdb.pl nt
```

### Compilation and Configuration
To install Pamir, you need to first fetch Pamir from our [git repository](https://github.com/vpc-ccg/pamir) or download the corresponding compressed files. 
```
git clone --recursive https://bitbucket.org/compbio/pamir.git
cd pamir
```

Then you need to update *pamir.config* in **pamir** folder with your paths for the binaries *samtools*, *velveth*, *velvetg*, *blastn* and the blast database folder. Open *pamir.config* using your preferred text editor, and modify your paths for the binaries:

```
SAMTOOLS=/your/path/to/samtools-1.3.1/samtools
VELVETH=/your/path/to/velvet/velveth
VELVETG=/your/path/to/velvet/velvetg
BLASTN=/your/path/to/ncbi-blast-2.5.0+/bin/blastn
BLASTDB=/your/path/to/ncbi-blast-2.5.0+/db/
```

Make the project
```
pamir$ make -j
```

Now you are ready to go!

### Troubleshooting
*SSE4 for mrsFAST* If sse4 is not supported in your system, you need to disable the flag of mrsfast by either modifying line 5 of pamir/Makefile to be
```
make with-sse4=no -C ../mrsfast
```

## Commands Options
Pamir is designed to detect novel sequence insertions based on one-end anchors (OEA) and orphans from paired-end Whole Genome Sequencing (WGS) reads.

Note that **reference genome** is required for running Pamir in addition to **sequencing or mapping** data.


A typical command to start Pamir is 
```
./pamir.py -p my_project -r ref.fa --files [filetype]=[filenames]
```

### Project Name
To run Pamir you have to specify a project name such that Pamir will create a folder to store the results and intermediate files. You need to specify project name by `-p`. 

### Data Preparation
#### Required Information
Two information are required for running Pamir:
1. Reference Genome: You need to provide the reference genome in single fasta file by specifying the parameters `-r` or `--reference`.

2. Masking File: You can provide a file for masking reference genome. For example,  you can ask Pamir to ignore events in repeat regions by giving ` -m repeat.mask` . When you only want to consider events in genic regions, use   `-m genic.region --invert-masker` and Pamir will mask those regions not in the given file.

#### Read Length
Now Pamir only accepts WGS datasets in which two mates of all reads are of **equal length**.

### Sequencing Data
Pamir can take either FASTQ and SAM files as its input. It has three different options to accept inputs:
1. **SAM/by mrsFAST-best-search**: A paired-end mapping result of your WGS data which satisfies the following conditions: 
   1. Two mates from a read are grouped together.
   2. All mates are of equal length.
    
   For example, a *best-mapping* SAM file is a valid input file for Pamir through ```--files mrsfast-best-search=wgs.sam```.

   You can also give multiple best-mapping files separated by comma, ```--files mrsfast-best-search=sample1.sam,sample2.sam,sample3.sam```,
   or the path to the folder that includes the inputs files, ```--files mrsfast-best-search=/directory/to/sample_best_mapping_sam_files/```.

2. **FASTQ**: Pamir also accepts FASTQ format as the input data once it is a single gzipped file in which two (equal-length) mates of a read locate consecutively. You can specify by giving `--files fastq=wgs.fastq.gz`.

3. **Alignment file SAM/BAM**: Pamir also accepts any other alignment output sorted by readname. Alignment output can be in SAM or BAM format. You can specify by `--files alignment=wgs.sam or --files alignment=wgs.bam`.



### mrsFAST Parameters
Pamir uses mrsFAST for multi-mapping the orphan and OEA reads obtained from the best-mapping output. You can give your own mrsFAST parameters or Pamir will use the default values. Some of the parameters you may want to update are :
1. **--mrsfast-n**: Maximum number of mapping loci of anchor of an OEA. Anchor with higher mapping location will be ignored. 0 for considering all mapping locations. 
(Default = 50)
2. **--mrsfast-threads**: Number of the threads used by mrsFAST-Ultra for mapping. (Default = 1)
3. **--mrsfast-errors**: Number of the errors allowed by mrsFAST-Ultra for mapping. In default mode Pamir does not give any error number to mrsFAST-Ultra, in which case it calculates the error value as 0.06 x readlength. (Default = -1; 0.06% of the read length)
4. **--mrsfast-index-ws**: Window size used by mrsFAST-Ultra for indexing the reference genome. (Default = 12)


### Other Parameters
1. **--num-worker**: Number of independent prediction jobs to be created. You can define this parameter according to your core number. (Default = 1)
2. **--resume**: Restart pipeline of an existing project from the stage that has not been completed yet.
3. **--assembler**: The assembler to be used in orphan assembly stage (Options: velvet, minia, sga. Default = velvet).
## Output Formats
Pamir generates a [VCF file](https://samtools.github.io/hts-specs/VCFv4.2.pdf) for detected novel sequence insertions. You can run genotyping for each sample after obtaining the VCF file by:

``` 
python genotyping.py projectFolder/insertions.out_wodups_filtered_setcov_PASS.sorted reference.fa.masked sample1_FASTQ_1.fq sample1_FASTQ_2.fq [readlength] [SAMPLENAME] [mrsfast-min] [mrsfast-max] [projectFolder] [TEMPLATE_LEN]
```

The default template length is set to be 1000.
## Example Commands
1. To start a new analysis from a mrsfast-best mapping result SAM file:
```
$ ./pamir.py -p my_project -r ref.fa --files mrsfast-best-search=sample.sam
```
2. make a pooled-run with multiple samples separated by comma: 
```
$ ./pamir.py -p my_project -r ref.fa  --files mrsfast-best-search=sample.sam,sample2.sam,sample3.sam
```
3. To make a pooled-run with multiple samples which are in a folder called SAMPLEFOLDER: 
```
$ ./pamir.py -p my_project -r ref.fa  --files mrsfast-best-search=SAMPLEFOLDER
```
4. To start from another mapping tool's alignment result SAM/BAM file:
```
$ ./pamir.py -p my_project -r ref.fa --files alignment=sample.bam
```
5. To start from a gzipped fastq file,
```
$ ./pamir.py -p my_project -r ref.fa --files fastq=sample.fastq.gz
```
6. To ignore regions in a mask file (e.g., repeat regions),
```
$ ./pamir.py -p my_project -r ref.fa -m repeat.txt --files mrsfast-best-search=sample.sam
```
7. To analyze events only in some regions of the reference genome (e.g., genic regions), 
```
$ ./pamir.py -p my_project -r ref.fa -m genic.region --invert-mask  --files mrsfast-best-search=sample.sam
```
8. To make sure that mrsFAST will not report the mapping locations of an OEA read more after the 30th location:
```
$ ./mistrvar.py -p my_project -r ref.fa --mrsfast-n 30 --files mrsfast-best-search=sample.sam
```
9. To specify the core number for mrsFAST during multi-mapping of OEAs: 
```
$ ./pamir.py -p my_project -r ref.fa  --mrsfast-threads 8 --files mrsfast-best-search=sample.sam
```
10. To speed up the prediction process by defining the independent prediction jobs according to available core numbers: 
```
$ ./pamir.py -p my_project -r ref.fa  --num-worker 20 --files mrsfast-best-search=sample.sam
```
11. To specify the assembler as sga for orphan assembly and also the number of prediction jobs will be 20: 
```
$ ./pamir.py -p my_project -r ref.fa  --num-worker 20 --assembler sga --files mrsfast-best-search=sample.sam
```
12. To resume from the previously finished stage: 
```
$ ./pamir.py -p my_project --resume
```

### Some Invalid Commands
The following parameters are not accepted by Pamir:
1. Project name is missing:
```
$./pamir.py -r ref.fa --files alignment=sample.sam
```
2. Reference genome file is missing:
```
$./pamir.py -p my_project --files alignment=sample.sam
```
3. Incorrect path of the mask file:
```
$./pamir.py -p my_project -m non-exist-mask-file --files alignment=sample.sam
```
4. No input sequencing files:
```
$./pamir.py -p my_project -r ref.fa
```
5. Multiple sequencing sources:
```
$./pamir.py -p my_project  -r ref.fa --files mrsfast-best-search=sample.sam fastq=sample2.fastq.gz
```

## Contact & Support

Feel free to drop any inquiry to [pinarkavak at gmail dot com]()    .