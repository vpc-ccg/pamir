#!/usr/bin/env python
import os, sys, errno, argparse, subprocess, fnmatch, ConfigParser, shutil

def usage():
	print '\nUsage: python genotyping.py VCF REF SEQ1.fastq SEQ2.fastq readlength EXTENSION mrsFAST-min mrsFAST-max workdir TLEN'
	sys.exit(-1)
def main():
	args = sys.argv[1:]
	if len(args) !=10:
		usage()
	REF=sys.argv[2]
	#FILE needs to be sorted according to 1st and 2nd columns
	FILE=sys.argv[1]
	SEQ1=sys.argv[3]
	SEQ2=sys.argv[4]
	readlen=int(sys.argv[5])
	EXT=sys.argv[6]
	MIN=sys.argv[7]
	MAX=sys.argv[8]
	DIR=sys.argv[9]
	TLEN = int(sys.argv[10])
	folder  ="{0}/genotype".format(DIR)
	os.system("mkdir -p {0}".format(folder))
	start=1
	start2=1
	f = open("{0}/allinsertions.fa".format(folder),"w")
	coor = open("{0}/allinsertions.coor".format(folder),"w")
	f.write(">1\n")
	f2 = open("{0}/allref.fa".format(folder),"w")
	coor2 = open("{0}/allref.coor".format(folder),"w")
	f2.write(">1\n")
	prev_loc=""
	prev_chrN=""
	failed = 0
	passed = 0
	a = 2
	vcfcontent = dict()
	fil = open (FILE + "_genotype_"+EXT,"w")
	os.system("samtools faidx {0}".format(REF))
	with open(FILE) as insertions:
		for line in insertions:
			elem_ins=line.split()
			if len(elem_ins) > 0:
				chrN    = elem_ins[0]
				loc     = elem_ins[1]
				length  = elem_ins[3]
				seq     = elem_ins[2]
				be=int(loc)-1-TLEN
				if be<0:
					be =0
				en=int(loc)+TLEN
				open("{0}/left.bed".format(folder),"w").write("{0}\t{1}\t{2}\n".format(chrN,be,int(loc)-1))
				open("{0}/right.bed".format(folder),"w").write("{0}\t{1}\t{2}\n".format(chrN,int(loc)-1,en))
				os.system("bedtools getfasta -bed {0}/left.bed -fi {1} -fo {0}/left.fa".format(folder, REF))
				os.system("bedtools getfasta -bed {0}/right.bed -fi {1} -fo {0}/right.fa".format(folder, REF))
				left=open("{0}/left.fa".format(folder),"r").readlines()[-1]
				right=open("{0}/right.fa".format(folder),"r").readlines()[-1]
				seqfaref = "{0}{1}".format(left[:len(left)-1], right[:len(right)-1] )
				seqfains = "{0}{1}{2}".format(left[:len(left)-1], seq, right[:len(right)-1] )
				f.write(seqfains)
				f2.write(seqfaref)
				loc2 =loc
				if(chrN==prev_chrN and loc==prev_loc):
					loc2 = loc+"_"+str(a)
					key = chrN + "_" + loc2
					vcfcontent[key]=[]
					vcfcontent[key].append(length)
					vcfcontent[key].append(seq)
					vcfcontent[key].append(0)
					vcfcontent[key].append(0)
					a+=1
				else:
					prev_loc =loc
					prev_chrN =chrN
					key = chrN + "_" + loc
					vcfcontent[key]=[]
					vcfcontent[key].append(length)
					vcfcontent[key].append(seq)
					vcfcontent[key].append(0)
					vcfcontent[key].append(0)
					a=2
				coor.write("{0}_{1}\t{2}\n".format(chrN, loc2, start ))
				coor2.write("{0}_{1}\t{2}\n".format(chrN, loc2, start2 ))
				start=start+len(left)-1+len(seq)+len(right)-1
				start2=start2+len(left)-1+len(right)-1
	f.close()
	f2.close()
	coor.close()
	coor2.close()
	os.system("./mrsfast --index {0}/allref.fa > {0}/mrsfast.index.log".format(folder))
	os.system("./mrsfast --search {0}/allref.fa --pe --min {4} --max {5} -n 50 --threads 24 -o {0}/seq.mrsfast.ref.{1}.sam -e 3 --seq1 {2} --seq2 {3} --disable-sam-header --disable-nohits > {0}/seq.mrsfast.ref.{1}.sam.log".format(folder, EXT, SEQ1, SEQ2, MIN, MAX))
	os.system("./recalibrate {0}/allref.coor {0}/seq.mrsfast.ref.{1}.sam {0}/seq.mrsfast.ref.{1}.recal.sam".format(folder,EXT))
	os.system("sort -k 3,3 -k 4,4n {0}/seq.mrsfast.ref.{1}.recal.sam > {0}/seq.mrsfast.ref.{1}.recal.sam.sorted".format(folder,EXT))
	msamlist = open("{0}/seq.mrsfast.ref.{1}.recal.sam.sorted".format(folder,EXT),"r")

	os.system("./mrsfast --index {0}/allinsertions.fa > {0}/mrsfast.index2.log".format(folder))
	os.system("./mrsfast --search {0}/allinsertions.fa --pe --min {4} --max {5} -n 50 --threads 24 -o {0}/seq.mrsfast.ins.{1}.sam -e 3 --seq1 {2} --seq2 {3} --disable-sam-header --disable-nohits > {0}/seq.mrsfast.ins.{1}.sam.log".format(folder, EXT, SEQ1, SEQ2, MIN, MAX))
	os.system("./recalibrate {0}/allinsertions.coor {0}/seq.mrsfast.ins.{1}.sam {0}/seq.mrsfast.ins.{1}.recal.sam".format(folder,EXT))
	os.system("sort -k 3,3 -k 4,4n {0}/seq.mrsfast.ins.{1}.recal.sam > {0}/seq.mrsfast.ins.{1}.recal.sam.sorted".format(folder,EXT))
	msamlist2 = open("{0}/seq.mrsfast.ins.{1}.recal.sam.sorted".format(folder, EXT),"r")
	i=0
	j=0
	chrName=""
	passNum =0 
	num =0
	breakpoint =TLEN+1
	last =""
	refsupport=0
	altsupport=0
	line = msamlist.readline()
	while(line !=''):
		refsupport=0
		ispass=0
		splitline = line.split()
		flag = int(splitline[1])
		locName = splitline[2]
		first_sep = locName.find("_")
		last_sep = locName.rfind("_")
		chrName = locName[0:first_sep]
		if(first_sep ==last_sep):
			last_sep = len(locName)
		firstloc = int(splitline[3])
		location = locName[first_sep+1:last_sep]
		if( flag & 2 == 2 and firstloc < breakpoint-10 and firstloc + readlen >= breakpoint+10):
			refsupport+=1
		i +=1;
		line = msamlist.readline()
		while(line !=''):
			splitline = line.split()
			flag = int(splitline[1])
			nextlocName = splitline[2]
			tmp = int(splitline[3])
			if (flag & 2 == 2 and tmp < breakpoint-10 and tmp + readlen >= breakpoint+10):
				refsupport+=1
			if(nextlocName != locName):
				break;
			lastloc = tmp
			i+=1
			line = msamlist.readline()
		vcfcontent[locName][2]=refsupport
	line2 = msamlist2.readline()
	while (line2 !=''):
		altsupport_left=0
		altsupport_right=0
		splitline2 = line2.split()
		flag2 = int(splitline[1])
		locName2 = splitline2[2]
		secondbreakpoint = breakpoint+int(vcfcontent[locName2][0])
		first_sep2 = locName2.find("_")
		last_sep2 = locName2.rfind("_")
		chrName2 = locName2[0:first_sep2]
		if(first_sep2 ==last_sep2):
			last_sep2 = len(locName2)
		firstloc2 = int(splitline[3])
		location2 = locName2[first_sep2+1:last_sep2]
		if(flag2 & 2 ==2 and firstloc2 < breakpoint-10 and firstloc2 + readlen >= breakpoint+10):
			altsupport_left+=1
		if(flag2 & 2 == 2 and firstloc2 < secondbreakpoint-10 and firstloc2 + readlen >= secondbreakpoint+10):
			altsupport_right+=1
		j +=1;
		line2 = msamlist2.readline()
		while(line2 != ''):
			splitline2 = line2.split()
			flag2 = int(splitline2[1])
			nextlocName2 = splitline2[2]
			tmp2 = int(splitline2[3])
			if ( flag2 & 2 == 2 and tmp2 < breakpoint-10 and tmp2 + readlen >= breakpoint+10):
				altsupport_left+=1
			if (flag2 & 2 == 2 and tmp2 < secondbreakpoint-10 and tmp2 + readlen >= secondbreakpoint+10):
				altsupport_right+=1
			if(nextlocName2 != locName2):
				break;
			lastloc2 = tmp2
			j+=1
			line2 = msamlist2.readline()
		altsupport = float((altsupport_left+altsupport_right))/2
		vcfcontent[locName2][3]= altsupport
	for a in vcfcontent:
		elem = vcfcontent[a]
		chrNameLoc=a.split('_')
		ratio=0
		if elem[2]==0 and elem[3]==0:
			fil.write("{0}\t{1}\t{2}\t{3}\t1/1\t{4}\t{5}\t{6}\n".format(chrNameLoc[0], chrNameLoc[1], elem[1], elem[0], elem[2], elem[3], str(ratio)) )
		else:
			ratio = (float(elem[3])-float(elem[2]))/(float(elem[3])+float(elem[2]))
			if(ratio >= 0.3):
				fil.write("{0}\t{1}\t{2}\t{3}\t1/1\t{4}\t{5}\t{6}\n".format(chrNameLoc[0], chrNameLoc[1], elem[1], elem[0], elem[2], elem[3], str(ratio)) )
			elif(ratio <=-0.3):
				fil.write("{0}\t{1}\t{2}\t{3}\t0/0\t{4}\t{5}\t{6}\n".format(chrNameLoc[0], chrNameLoc[1], elem[1], elem[0], elem[2], elem[3], str(ratio)) )
			else:
				fil.write("{0}\t{1}\t{2}\t{3}\t0/1\t{4}\t{5}\t{6}\n".format(chrNameLoc[0], chrNameLoc[1], elem[1], elem[0], elem[2], elem[3], str(ratio)) )
	fil.close()

#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

