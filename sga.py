#!/usr/bin/env python
import os, sys, errno, argparse, subprocess, fnmatch, ConfigParser, datetime

def main():
	SGA="/cs/compbio3/yenyil/Pinar/LIB/sgabin/bin/sga"
	FASTQ  = os.path.abspath(sys.argv[1]);
	print "sga start " + FASTQ + " " + str(datetime.datetime.now()) 
	dirname = FASTQ.split("/");
	diir= ""
	for i in range(0,len(dirname)-1):
		diir+=dirname[i]+"/"
	os.chdir(diir)
	os.system(SGA + " preprocess " + FASTQ + " -o " + FASTQ + ".preprocess");
	os.system(SGA + " index -a ropebwt -t 8 " + FASTQ + ".preprocess");
	os.system(SGA + " correct -a overlap -t 8 " + FASTQ + ".preprocess");
	os.system(SGA + " index -a ropebwt -t 8 " + FASTQ + ".ec.fa");
	os.system(SGA + " filter -t 8 " + FASTQ + ".ec.fa --no-kmer-check");
	os.system(SGA + " overlap -m 70 " + FASTQ + ".ec.filter.pass.fa -t 8");
	os.system(SGA + " assemble -o " + FASTQ + " -m 70 " + FASTQ + ".ec.filter.pass.asqg.gz");
	fin = open(FASTQ + "-contigs.fa","r")
	fout = open(FASTQ + ".sgaout.fa","w")
	contigs = fin.readlines()
	i =0 
	while i < len(contigs):
		splitted=contigs[i].split(" ")
		if int(splitted[1]) > 100:
			fout.write(splitted[0]+ "\n" +contigs[i+1])
		i+=2
	fin.close()
	fout.close()

#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

