#!/usr/bin/env python
#-*- encoding:utf-8 -*-
import sys
import math
import os
class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def usage():
	print '\nUsage: python extractContaminants.py contaminationFile orphan.fq.contigs.fa orphan.fq.contigs.wocontamination.fa'
	sys.exit(-1)

def main():
	args = sys.argv[1:]
	if len(args) !=3:
		usage()

	fcontamination = open(sys.argv[1], 'r')
	fcontig = open(sys.argv[2], 'r')
	fout = open(sys.argv[3],'w')
	
	contamination = fcontamination.readline()
	while contamination!='':
		contamination = contamination.split()[0]
		contig = fcontig.readline()
		if(contig!=''):
			contig = contig.strip().split(">")[1]
			while contamination !=contig and contig !='':
				fout.write(">" + contig + "\n")
				fout.write(fcontig.readline())
				contig = fcontig.readline()
				if(contig!=''):
					contig = contig.strip().split(">")[1]
			if(contamination==contig and contig!=''):
				fcontig.readline()
			contamination = fcontamination.readline()
	contig = fcontig.readline()
	while contig !='':
		fout.write(contig)
		fout.write(fcontig.readline())
		contig = fcontig.readline()
	fout.close()

if __name__ == "__main__":
	sys.exit(main())
