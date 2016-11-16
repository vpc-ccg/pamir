#!/usr/bin/env python
#-*- encoding:utf-8 -*-
import sys
import math
import os
class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg
def usage():
	print '\nUsage: python generateClustersForSetCover.py VCF SAM forSetCover'
	sys.exit(-1)

def main():
	args = sys.argv[1:]
	if len(args) !=3:
		usage()

	vcffile = open(sys.argv[1],"r").readlines()
	samfile = open(sys.argv[2],"r").readlines()
	fout = open(sys.argv[3],"w")
	
	vcfnum = 0
	samline = 0
	readnames = set()
	clusterId =1
	samel = samfile[samline].strip().split()
	samline+=1
	while vcfnum < len(vcffile):
		readnames.clear()
		elem = vcffile[vcfnum].strip().split()
		print elem[0]
		while vcfnum< len(vcffile) and elem[5]=="FAIL":
			vcfnum+=1
			if vcfnum < len(vcffile):
				elem = vcffile[vcfnum].strip().split()
			else:
				break
		name = elem[0]
		if(elem[0].count("-")>1):
			name = elem[0][0:elem[0].rfind("-")]

		while vcfnum < len(vcffile) and samline < len(samfile) and elem[0]>samel[2]:
			samel = samfile[samline].strip().split()
			samline+=1
		if elem[0]==samel[2]:
			while vcfnum < len(vcffile) and samline < len(samfile) and elem[0]==samel[2]:
				readnames.add(samel[0])
				samel = samfile[samline].strip().split()
				samline+=1

		if vcfnum+1 < len(vcffile):
			nextelem = vcffile[vcfnum+1].strip().split()
		print nextelem[0]
		while vcfnum+1 < len(vcffile) and nextelem[0].count("-") >1 and nextelem[0][0:nextelem[0].rfind("-")]==name and nextelem[5]=="FAIL":
			vcfnum+=1
			if(vcfnum+1< len(vcffile)):
				nextelem = vcffile[vcfnum+1].strip().split()
		print nextelem[0]
		while vcfnum< len(vcffile) and nextelem[0].count("-") >1 and nextelem[0][0:nextelem[0].rfind("-")]==name:
			while vcfnum < len(vcffile) and samline < len(samfile) and nextelem[0]==samel[2]:
				readnames.add(samel[0])
				samel = samfile[samline].strip().split()
				samline+=1
			vcfnum+=1
			if(vcfnum < len(vcffile)):
				nextelem = vcffile[vcfnum+1].strip().split()
			while(vcfnum < len(vcffile) and nextelem[5]=="FAIL"):
				vcfnum+=1
				if(vcfnum+1 < len(vcffile)):
					nextelem = vcffile[vcfnum+1].strip().split()

		if len(readnames)>0:
			fout.write(str(clusterId)+ " "+ str(len(readnames))+ " "+name+"\n")
			clusterId+=1
			for a in readnames:
				fout.write(a+"\n")
		vcfnum+=1

	fout.close()

if __name__ == "__main__":
	sys.exit(main())
