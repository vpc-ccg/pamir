#!/usr/bin/env python
#-*- encoding:utf-8 -*-
import sys
import math
import os
class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def usage():
	print('\nUsage: python extractContaminants.py contaminationFile orphan.fq.contigs.fa orphan.fq.contigs.wocontamination.fa')
	sys.exit(-1)

def main():
	args = sys.argv[1:]
	if len(args) !=3:
		usage()

	cont = {}
	sr = open(sys.argv[1], 'r')
	for line in sr:
		cont_id = line.strip().split()[0]
		cont[cont_id] = 1
	sr.close()



	fcontig = open(sys.argv[2], 'r')
	fout = open(sys.argv[3],'w')
	
	sw = open( sys.argv[3], 'w')
	sr = open( sys.argv[2], 'r')
	flag = 1
	for line in sr:
		if ">" == line[0]:
			flag = 1
			c_id = line.strip().split()[0][1:]
			if c_id in cont:
				flag = 0
		if 1 == flag:
			sw.write(line)
	sw.close()
	sw.close()

if __name__ == "__main__":
	sys.exit(main())
