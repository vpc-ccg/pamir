#!/usr/bin/env python
import os, sys, errno, argparse, subprocess, fnmatch, ConfigParser, shutil
def main():
	REF=sys.argv[2]
	FILE=sys.argv[1]
	SEQ1=sys.argv[3]
	SEQ2=sys.argv[4]
	readlen=int(sys.argv[5])
	if(os.path.isfile(FILE+"_filter_genotype_recal")):
		os.unlink(FILE+"_filter_genotype_recal")
	folder  ="filter_genotype"
	os.system("mkdir -p {0}".format(folder))
	start=1
	start2=1
	f = open("{0}/allinsertions.fa".format(folder),"w")
	coor = open("{0}/allinsertions.coor".format(folder),"w")
	f2 = open("{0}/allref.fa".format(folder),"w")
	coor2 = open("{0}/allref.coor".format(folder),"w")
	f.write(">1\n")
	f2.write(">1\n")
	prev_loc=""
	failed = 0
	passed = 0
	a = 2
	vcfcontent = dict()
	fil = open (FILE + "_filter_genotype_recal","w")
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
				rightCl = int(1000+int(length))
				be=int(loc)-1-1000
				if be < 0:
					be = 0
				en=int(loc)+1000
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
				if(loc==prev_loc):
					loc2 = loc+"_"+str(a)
					key = chrN + "_" + loc2
					vcfcontent[key]=[]
					vcfcontent[key].append(length)
					vcfcontent[key].append(seq)
					vcfcontent[key].append(seqfains)
					a+=1
				else:
					prev_loc =loc
					key = chrN + "_" + loc
					vcfcontent[key]=[]
					vcfcontent[key].append(length)
					vcfcontent[key].append(seq)
					vcfcontent[key].append(seqfains)
					a=2
				coor.write("{0}_{1}\t{2}\n".format(chrN, loc2, start ))
				coor2.write("{0}_{1}\t{2}\n".format(chrN, loc2, start2 ))
				start=start+len(left)-1+len(seq)+len(right)-1
				start2=start2+len(left)-1+len(right)-1
	f.close()
	f2.close()
	coor.close()
	coor2.close()
	
	os.system("../mrsfast --index {0}/allinsertions.fa > {0}/mrsfast.index2.log".format(folder))
	os.system("../mrsfast --search {0}/allinsertions.fa --pe --min 300 --max 600 -o {0}/seq.mrsfast.ins.sam -e 3 --seq1 {1} --seq2 {2} --disable-sam-header --disable-nohits > {0}/seq.mrsfast.ins.sam.log".format(folder, SEQ1, SEQ2))
	os.system("../recalibrate {0}/allinsertions.coor {0}/seq.mrsfast.ins.sam {0}/seq.mrsfast.ins.recal.sam".format(folder))
	os.system("sort -k 3,3d -k 4,4n {0}/seq.mrsfast.ins.recal.sam > {0}/seq.mrsfast.ins.recal.sam.sorted".format(folder))
	msamlist = open("{0}/seq.mrsfast.ins.recal.sam.sorted".format(folder),"r").readlines()

	os.system("../mrsfast --index {0}/allref.fa > {0}/mrsfast.index.log".format(folder))
	os.system("../mrsfast --search {0}/allref.fa --pe --min 300 --max 600 -o {0}/seq.mrsfast.ref.sam -e 3 --seq1 {1} --seq2 {2} --disable-sam-header --disable-nohits > {0}/seq.mrsfast.ref.sam.log".format(folder, SEQ1, SEQ2))
	os.system("../recalibrate {0}/allref.coor {0}/seq.mrsfast.ref.sam {0}/seq.mrsfast.ref.recal.sam".format(folder))
	os.system("sort -k 3,3d -k 4,4n {0}/seq.mrsfast.ref.recal.sam > {0}/seq.mrsfast.ref.recal.sam.sorted".format(folder))
	msamlist2 = open("{0}/seq.mrsfast.ref.recal.sam.sorted".format(folder),"r").readlines()

	i=0
	j=0
	chrName=""
	num =0
	breakpoint =1001
	while(i < len(msamlist) and j < len(msamlist2)):
		lsupport = 0
		rsupport = 0
		tsupport = 0
		refsupport=0
		altsupport=0
		msamsplit = msamlist[i].split()
		locName = msamsplit[2]
		first_sep = locName.find("_")
		last_sep = locName.rfind("_")
		chrName = locName[0:first_sep]
		if(first_sep ==last_sep):
			last_sep = len(locName)
		location = locName[first_sep+1:last_sep]
		firstloc = int(msamsplit[3])
		rightCl  = int(1000+int(vcfcontent[locName][0])+1)
		if firstloc <= leftCl:
			lsupport += 1
		if firstloc >= rightCl:
			rsupport += 1
		if(firstloc < breakpoint-10 and firstloc + readlen >= breakpoint+10):
			altsupport+=1
		i +=1
		while (i < len(msamlist)):
			msamsplit = msamlist[i].split()
			nextlocName = msamsplit[2]
			tmp = int (msamsplit[3])
			if (tmp < breakpoint-10 and tmp + readlen >= breakpoint+10):
				altsupport+=1
			if(nextlocName != locName):
				end = i
				break;
			lastloc = tmp
			if(lastloc <= leftCl):
				lsupport+=1
			if(lastloc >= rightCl):
				rsupport+=1
			i+=1
		tsupport=lsupport+rsupport
		ispass=0
		if(firstloc<=leftCl and lastloc >=rightCl):
			ispass=1
		else:
			reffastaleft = vcfcontent[locName][2][0:leftCl]
			if(reffastaleft.count('N') >= (len(reffastaleft)/2)):
				if(int(firstloc) <= leftCl or int(lastloc) >= rightCl):
					ispass =1
			reffastaright = vcfcontent[locName][2][rightCl:len(vcfcontent[locName][2])]
			if(reffastaright.count('N')>=(len(reffastaright)/2)):
				if(int(firstloc)<=leftCl or int(lastloc) >= rightCl):
					ispass =1
		if(ispass==1):
			msamsplit2 = msamlist2[j].split()
			locName2 = msamsplit2[2]
			if locName2==locName:
				location2 = location
				firstloc2 = int(msamsplit2[3])
				print firstloc2
				if (firstloc2 <= breakpoint - 10 and firstloc2 + readlen >= breakpoint + 10):
					refsupport+=1
				j+=1
				while j < len(msamlist2):
					msamsplit2 = msamlist2[j].split()
					nextlocName2 = msamsplit2[2]
					tmp2 = int(msamsplit2[3])
					print tmp2
					if(tmp2 < breakpoint-10 and tmp2 + readlen >= breakpoint +10):
						refsupport +=1
					if nextlocName2 != locName2: 
						end2 = j
						break;
					lastloc2 = tmp2
					j+=1

				print refsupport, altsupport
				if refsupport !=0 or altsupport != 0:
					elem = vcfcontent[locName]
					ratio = float(altsupport-refsupport)/float(altsupport+refsupport)
					if(ratio >= 0.3):
						fil.write("{0}\t{1}\t{2}\t{3}\t1/1\t{4}\t{5}\t{6}\n".format(chrName, location, elem[1], elem[0], refsupport, altsupport, str(ratio)))
					elif (ratio <= -0.3):
						fil.write("{0}\t{1}\t{2}\t{3}\t0/0\t{4}\t{5}\t{6}\n".format(chrName, location, elem[1], elem[0], refsupport, altsupport, str(ratio)))
					else:
						fil.write("{0}\t{1}\t{2}\t{3}\t0/1\t{4}\t{5}\t{6}\n".format(chrName, location, elem[1], elem[0], refsupport, altsupport, str(ratio)))
				elif refsupport==0 and altsupport ==0 :
					fil.write("{0}\t{1}\t{2}\t{3}\t0/00\t{4}\t{5}\t{6}\n".format(chrName, location, elem[1], elem[0], refsupport, altsupport, str(ratio)))
					
	fil.close()
						
	#os.system("grep PASS sniper_part_updated.vcf_sorted_wodups_filtered | awk '{print $2\"\t\"$4;}' | sort -k 1,1n > sniper_part_updated.vcf_sorted_wodups_filtered_PASS_loc")

#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

