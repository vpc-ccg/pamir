#Pamir: Insertion Discovery Tool for Whole Genome Sequencing Data

---

Pamir is a computational tool for detecting novel sequence insertions in single or multiple paired-end WGS(Whole Genome Sequencing) Illumina reads based on orphans and one-end anchors (OEA).

### Prerequisities and Compilation

Pamir relies on specific version of the following tools:  

g++ 4.9.0 or higher  
(https://gcc.gnu.org/releases.html)  

Python 2.7 or higher (needed for the package argparse )  

boost library 1.62 or higher  
(https://sourceforge.net/projects/boost/?source=directory)  
You also need to define the boost library by typing on your shell  
export BOOST INCLUDE= the/BOOST/version/include/ (directory of BOOST in your machine).

velvet 1.2.10 or higher

BLAST 2.3.0+ or higher

Latest BLAST nt database is also needed to be downloaded in dir/to/blast/db (needed for contamination detection).  
mkdir dir/to/blast/db  
cd dir/to/blast/db  
../bin/update blastdb.pl nt  

Then, just clone our repository:

```
git clone --recursive git@bitbucket.org:compbio/pamir.git
cd pamir
```

You need to link the velvet executables and blast directory inside pamir folder:  
ln -s /dir/to/velvet/velveth /dir/to/pamir/velveth  
ln -s /dir/to/velvet/velvetg /dir/to/pamir/velvetg  
ln -s /dir/to/blast folder /dir/to/pamir/blast  

Then make
```
pamir # make
```

### How do I run Pamir?
You can use 
```
pamir/pamir.py -h
```
to get a description of each parameter. For more details, please check doc/pamir_manual.pdf.


---


### Contact & Support

Feel free to drop any inquiry to [pinarkavak at gmail dot com](mailto:).
