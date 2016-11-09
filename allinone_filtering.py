#!/usr/bin/env python
import os, sys, errno, argparse, subprocess, fnmatch, ConfigParser, shutil
def main():
	REF			=	sys.argv[2]
	FILE		=	sys.argv[1]
	readlength  =	int(sys.argv[3])
	MIN			=	sys.argv[4]
	MAX			=	sys.argv[5]
	workdir		= 	sys.argv[6]
#	if(os.path.isfile(FILE+"_filtered_recal")):
#		os.unlink(FILE+"_filtered_recal")
	folder  ="{0}/filtering".format(workdir)
	os.system("mkdir -p {0}".format(folder))
	start=1
	f = open("{0}/allinsertions.fa".format(folder),"w")
	coor = open("{0}/allinsertions.coor".format(folder),"w")
	f.write(">1\n")
	prev_loc=""
	prev_chrN=""
	failed = 0
	passed = 0
	a = 2
	vcfcontent = dict()
	fil = open (FILE + "_filtered_recal","w")
	fil2 = open (FILE + "_filtered_recal_forSETCOVER","w")
	os.system("samtools faidx {0}".format(REF))
	with open(FILE) as insertions:
		for line in insertions:
			elem_ins=line.split()
			if len(elem_ins) > 0:
				chrN    = elem_ins[0]
				loc     = elem_ins[1]
				length  = elem_ins[2]
				seq     = elem_ins[3]
				leftCl  = int(900)
				rightCl = int(1000+int(length)+1)
				be=int(loc)-1-1001
				if be<0:
					be = 0
				en=int(loc)+1000
				open("{0}/left.bed".format(folder),"w").write("{0}\t{1}\t{2}\n".format(chrN,be,int(loc)-1))
				open("{0}/right.bed".format(folder),"w").write("{0}\t{1}\t{2}\n".format(chrN,int(loc)-1,en))
				os.system("bedtools getfasta -bed {0}/left.bed -fi {1} -fo {0}/left.fa".format(folder, REF))
				os.system("bedtools getfasta -bed {0}/right.bed -fi {1} -fo {0}/right.fa".format(folder, REF))
				left=open("{0}/left.fa".format(folder),"r").readlines()[-1]
				right=open("{0}/right.fa".format(folder),"r").readlines()[-1]
				seqfa = "{0}{1}{2}".format(left[:len(left)-1], seq, right[:len(right)-1] )
				f.write(seqfa)
				loc2 =loc
				if(chrN==prev_chrN and loc==prev_loc):
					loc2 = loc+"_"+str(a)
					key = chrN + "_" + loc2
					vcfcontent[key]=[]
					vcfcontent[key].append(length)
					vcfcontent[key].append(seq)
					vcfcontent[key].append(seqfa)
					a+=1
				else:
					prev_loc =loc
					prev_chrN =chrN
					key = chrN + "_" + loc
					vcfcontent[key]=[]
					vcfcontent[key].append(length)
					vcfcontent[key].append(seq)
					vcfcontent[key].append(seqfa)
					a=2
				coor.write("{0}_{1}\t{2}\n".format(chrN, loc2, start ))
				start=start+len(left)-1+len(seq)+len(right)-1
	f.close()
	coor.close()
	os.system("./mrsfast --index {0}/allinsertions.fa > {0}/mrsfast.index.log".format(folder))
	os.system("./mrsfast --search {0}/allinsertions.fa --pe --min {1} --max {2} -o {0}/seq.mrsfast.sam -e 3 --seq {3}/all_interleaved.fastq --threads 8 --disable-sam-header --disable-nohits > {0}/.seq.mrsfast.sam.log".format(folder, MIN, MAX, workdir))
	os.system("./recalibrate {0}/allinsertions.coor {0}/seq.mrsfast.sam {0}/seq.mrsfast.recal.sam".format(folder))
#	os.system("../sniper sort {0}/seq.mrsfast.recal.sam {0}/seq.mrsfast.recal.sam.sorted".format(folder))
	os.system("sort -k 3,3 -k 4,4n {0}/seq.mrsfast.recal.sam > {0}/seq.mrsfast.recal.sam.sorted".format(folder))
	msamlist = open("{0}/seq.mrsfast.recal.sam.sorted".format(folder),"r").readlines()
	i=0 
	chrName=""
	passNum =0 
	num =0
	while(i < len(msamlist)):
		lsupport=0
		rsupport=0
		tsupport=0
		splitmsam = msamlist[i].split()
		flag = int(splitmsam[1])
		locName = splitmsam[2]
		first_sep = locName.find("_")
		last_sep = locName.rfind("_")
		chrName = locName[0:first_sep]
		if(first_sep ==last_sep):
			last_sep = len(locName)
		location = locName[first_sep+1:last_sep]
		firstloc = int(splitmsam[3])
		tlen	 = int(splitmsam[8])
		rightCl = int(1000+int(vcfcontent[locName][0])+1)
		if flag & 2 ==2 and firstloc <=leftCl and firstloc + tlen >= 1001:
			lsupport+=1
		if flag & 2 == 2 and firstloc >=rightCl and firstloc + tlen + readlength <= rightCl:
			rsupport+=1
		i +=1;
		while(i < len(msamlist)):
			splitmsam	= msamlist[i].split()
			flag 		= int(splitmsam[1])
			nextlocName = splitmsam[2]
			tmp = int(splitmsam[3])
			tlen	 = int(splitmsam[8])
			if(nextlocName != locName):
				end = i
				break;
			lastloc = tmp
			if flag & 2 == 2 and lastloc <= leftCl and firstloc + tlen >= 1001:
				lsupport+=1
			if flag & 2 == 2 and lastloc >= rightCl and firstloc + tlen + readlength <= rightCl:
				rsupport+=1
			i+=1
		tsupport=lsupport+rsupport
		ispass 		   = 0	
		#if(firstloc<=leftCl and lastloc >= rightCl):
		if(lsupport > 0 and rsupport > 0):
			ispass =1
			passNum+=1
		else:
			reffastaleft=vcfcontent[locName][2][0:leftCl]
			if(reffastaleft.count('N')>=(len(reffastaleft)/2)):
				#if(int(firstloc)<=leftCl or int(lastloc)>=rightCl):
				if(lsupport > 0 or rsupport > 0):
					ispass=1
			reffastaright = vcfcontent[locName][2][rightCl:len(vcfcontent[locName][2])]
			if(reffastaright.count('N')>=(len(reffastaright)/2)):
				#if(int(firstloc)<=leftCl or int(lastloc)>=rightCl):
				if(lsupport > 0 or rsupport > 0):
					ispass=1
		num+=1
		elem = vcfcontent[locName]
		if ispass ==1:
			fil.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tPASS\n".format(chrName, location, elem[1], elem[0], lsupport, rsupport, tsupport) )
			fil2.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tPASS\n".format(locName, elem[1], elem[0], lsupport, rsupport, tsupport) )
			passed+=1
		else:
			fil.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tFAIL\n".format(chrName, location, elem[1], elem[0], lsupport, rsupport, tsupport) )
			fil2.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\tFAIL\n".format(locName, elem[1], elem[0], lsupport, rsupport, tsupport) )
			failed+=1
	fil.close()
	fil2.close()
	os.system("grep PASS "+FILE+"_filtered_recal | awk '{print $2\"\t\"$4;}' | sort -k 1,1n > "+FILE+"_filtered_recal_PASS_loc")
	os.system("perl sort_file.pl "+FILE+"_filtered_recal_forSETCOVER > "+FILE+"_filtered_recal_forSETCOVER.sorted")

#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

