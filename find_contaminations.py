#!/usr/bin/env python
#-*- encoding:utf-8 -*-
import sys
import math
import os
class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def usage():
	print '\nUsage: python contaminantFinder.py input output'
	sys.exit(-1)

def main():
	args = sys.argv[1:]
	if len(args) !=2:
		usage()

	fin = open(sys.argv[1], 'r')
	fout = open(sys.argv[2], 'w')
	contiglist=[]
	contiginfo=[]
	notcontaminant = ['Human', 'H.sapiens', 'Homo', 'sapiens', 'troglodytes', 'Gorilla', 'pygmaeus', 'mulatta', 'Rhinopithecus', 'Papio', 'Macaque', 'musculus', 'Callithrix', 'Colobus', 'fascicularis', 'Ailuropoda', 'Aotus', 'abelii', 'Pongo', 'Microcebus', 'Chinchilla', 'Macaca', 'macaque', 'aethiops', 'Mouse', 'canadensis','monkey','squirrel', 'Pig', 'Canadensis', 'Ovis', 'Chlorocebus','Myotis','gallus','Gallus','citb_109_p_11,', 'Methylobacterium']
	i=0
	line = fin.readline()
	while(line!=''):
		elems = line.strip().split();
		contiginfo=[]
		contiginfo.append(elems[1]);
		line = fin.readline()
		elems = line.strip().split("=");
		contiginfo.append(elems[1]);
		line = fin.readline()
		if line!='':
			prev_line = line;
			while( line !='' and not(line[0]=="Q" and line[1]=="u" and line[2]=="e" and line[3]=="r" and line[4]=="y")):
				elems= line.strip().split();
				contiginfo.append(elems[1]+" "+elems[2]);
				line = fin.readline()
			if(prev_line!= line):
				contiglist.append(contiginfo);
	for i in range(0,len(contiglist)):
		contigel = contiglist[i];
		f = 2;
		notprint = 0
		while f<15 and f < len(contigel):
			if(contigel[f].split()[0] in notcontaminant or contigel[f].split()[1] in notcontaminant):
				notprint = 1
			f+=1
		if(notprint==0):
			toprint=""
			for a in contigel:
				toprint+= a+"\t"
			fout.write("{0}\n".format(toprint)); 
		
	fout.close()

if __name__ == "__main__":
	sys.exit(main())
