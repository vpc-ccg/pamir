# **Pamir**: Insertion Discovery Tool for Whole Genome Sequencing Data

---

Pamir is a computational tool for detecting novel sequence insertions in single or multiple paired-end WGS(Whole Genome Sequencing) Illumina reads based on orphans and one-end anchors (OEA).

## Prerequisities and Compilation
---

Pamir relies on specific version of the following tools:  

g++ 4.9.0 or higher  
(https://gcc.gnu.org/releases.html)  

Python 2.7 or higher (needed for the package argparse )  

boost library 1.62 or higher  
(https://sourceforge.net/projects/boost/?source=directory) 

You also need to define the boost library by typing on your shell  
export BOOST\_INCLUDE= the/BOOST/version/include/ (directory of BOOST in your machine).

velvet 1.2.10 or higher

BLAST 2.3.0+ or higher

Latest BLAST nt database is also needed to be downloaded in dir/to/blast/db (needed for contamination detection).  
mkdir dir/to/blast/db  
cd dir/to/blast/db  
../bin/update blastdb.pl nt  

Then, just clone our repository:

```
git clone --recursive https://bitbucket.org/compbio/pamir.git
cd pamir
```

You need to update pamir.config in pamir folder with your paths for the binaries samtools, velveth, velvetg, blastn and the blast database folder db: 
```
vim pamir.config
Write your paths for the binaries:
SAMTOOLS=/your/path/to/samtools-1.3.1/samtools
VELVETH=/your/path/to/velvet/velveth
VELVETG=/your/path/to/velvet/velvetg
BLASTN=/your/path/to/ncbi-blast-2.5.0+/bin/blastn
BLASTDB=/your/path/to/ncbi-blast-2.5.0+/db/
```

Then make
```
pamir$ make
```

## How do I run Pamir?
---

You can use either of these commands:
```
pamir$ ./pamir.py -h
pamir$ python pamir.py -h
```
to get a description of each parameter. For more details, please check doc/pamir_manual.pdf.


---


## Contact & Support
---

Feel free to drop any inquiry to [pinarkavak at gmail dot com](mailto:).
