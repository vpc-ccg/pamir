#!/usr/bin/env python
import os, sys, errno, argparse, subprocess, fnmatch, ConfigParser, shutil


def usage():
	print '\nUsage: python filtering.py VCF REF mrsFAST-min mrsFAST-max workdir TLEN'
	sys.exit(-1)

##############################
def load_fasta( fasta_file ):
	ref_list = [] # keep order of reference in original input file
	ref_dict = {}
	tmp_list = []
	ref_id   = ""
	sr = open(fasta_file, "r")
	for line in sr:
		if ( ">" == line[0]):
			if ( ref_id != "" ):
				ref_dict[ ref_id] ="".join(tmp_list)
			ref_id   = line.strip().split()[0][1:]
			if ref_id not in ref_dict:
				ref_list.append( ref_id )
			del tmp_list[:]
		else:
			tmp_list.append( line.strip() )
	
	# adding records for the last chromosome
	if ( ref_id != "" ):
		ref_dict[ ref_id ] = "".join(tmp_list)

	sr.close()

	return ref_dict, ref_list
################################
def get_bed_seq( ref_dict, ref, start, end):
	if ref not in ref_dict:
		print "Error: " +ref + " not found in reference genome"
		exit(-1)
	return ref_dict[ref][start:end]
################################
def main():
	args = sys.argv[1:]
	if len(args) !=6:
		usage()
	MRSFAST= "mrsfast/mrsfast"
	REF			=	sys.argv[2]
	FILE		=	sys.argv[1]
	#mrsfast min
	MIN			=	sys.argv[3]
	#mrsfast max
	MAX			=	sys.argv[4]
	workdir		= 	sys.argv[5]
	#how many bp before the breakpoint and after the breakpoint on ref to get. (1000 for now)
	TLEN		=	int(sys.argv[6])
	tmpf = open("{0}/all_interleaved.fastq".format(workdir),"r")
	tmp = tmpf.readline()
	tmp = tmpf.readline().strip()
	readlength=int(len(tmp))
	tmpf.close()
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
	fil = open (FILE + "_filtered","w")
	fil2 = open (FILE + "_filtered_for_setcov","w")

	ref_dict, ref_list = load_fasta(REF)



	with open(FILE) as insertions:
		for line in insertions:
			elem_ins=line.split()
			if len(elem_ins) > 0:
				chrN    = elem_ins[0]
				loc     = elem_ins[1]
				r2		= elem_ins[2]
				r3		= elem_ins[3]
				r4 		= elem_ins[4]
				r5 		= elem_ins[5]
				r6 		= elem_ins[6]
				r7 		= elem_ins[7]
				tmplst 	= r7.split(";")
				length 	= int(tmplst[1].split("=")[1])
				seq		= tmplst[5].split("=")[1]
				leftCl  = int(TLEN-readlength)
				rightCl = int(TLEN+length+1)
				be=int(loc)-1-(TLEN+1)
				if be<0:
					be = 0
				en=int(loc)+TLEN
				left  = get_bed_seq( ref_dict, chrN, be, int(loc)-1)
				right = get_bed_seq( ref_dict, chrN, int(loc)-1, en)
				seqfa = left + seq + right
				f.write(seqfa)
				loc2 =loc
				if(chrN==prev_chrN and loc==prev_loc):
					loc2 = loc+"-"+str(a)
					key = chrN + "-" + loc2
					vcfcontent[key]=[]
					vcfcontent[key].append(length)
					vcfcontent[key].append(seq)
					vcfcontent[key].append(seqfa)
					vcfcontent[key].append(loc)
					vcfcontent[key].append(r2)
					vcfcontent[key].append(r3)
					vcfcontent[key].append(r4)
					vcfcontent[key].append(r5)
					vcfcontent[key].append(r6)
					vcfcontent[key].append(r7)
					a+=1
				else:
					prev_loc =loc
					prev_chrN =chrN
					key = chrN + "-" + loc
					vcfcontent[key]=[]
					vcfcontent[key].append(length)
					vcfcontent[key].append(seq)
					vcfcontent[key].append(seqfa)
					vcfcontent[key].append(loc)
					vcfcontent[key].append(r2)
					vcfcontent[key].append(r3)
					vcfcontent[key].append(r4)
					vcfcontent[key].append(r5)
					vcfcontent[key].append(r6)
					vcfcontent[key].append(r7)
					a=2
				coor.write("{0}-{1}\t{2}\n".format(chrN, loc2, start ))
				start=start+len(left)+len(seq)+len(right)
	f.close()
	coor.close()
	os.system(MRSFAST + " --index {0}/allinsertions.fa > {0}/mrsfast.index.log".format(folder))
	os.system(MRSFAST + " --search {0}/allinsertions.fa --pe --min {1} --max {2} -o {0}/seq.mrsfast.sam -e 3 --seq {3}/all_interleaved.fastq --threads 72 --disable-sam-header --disable-nohits > {0}/.seq.mrsfast.sam.log".format(folder, MIN, MAX, workdir))
	os.system("./recalibrate {0}/allinsertions.coor {0}/seq.mrsfast.sam {0}/seq.mrsfast.recal.sam".format(folder))
	os.system("sort -k 3,3 -k 4,4n {0}/seq.mrsfast.recal.sam > {0}/seq.mrsfast.recal.sam.sorted".format(folder))
	msamlist = open("{0}/seq.mrsfast.recal.sam.sorted".format(folder),"r")
	chrName=""
	passNum =0 
	line = msamlist.readline()
	while(line!=''):
		lsupport=0
		rsupport=0
		tsupport=0
		splitmsam = line.split()
		flag = int(splitmsam[1])
		locName = splitmsam[2]
		first_sep = locName.find("-")
		last_sep = locName.rfind("-")
		chrName = locName[0:first_sep]
		if(first_sep ==last_sep):
			last_sep = len(locName)
		location = locName[first_sep+1:last_sep]
		firstloc = int(splitmsam[3])
		tlen	 = int(splitmsam[8])
		rightCl = int(TLEN+int(vcfcontent[locName][0])+1)
		if flag & 2 == 2:
			if firstloc <=leftCl and firstloc + tlen >= (TLEN+1):
				lsupport+=1
			elif firstloc >=rightCl and firstloc + tlen + readlength <= rightCl:
				rsupport+=1
		line = msamlist.readline()
		while(line !=''):
			splitmsam	= line.split()
			flag 		= int(splitmsam[1])
			nextlocName = splitmsam[2]
			nextloc = int(splitmsam[3])
			tlen	 = int(splitmsam[8])
			if(nextlocName != locName):
				break;
			if flag & 2 ==2:
				if nextloc <= leftCl and nextloc + tlen >= (TLEN+1):
					lsupport+=1
				elif nextloc >= rightCl and nextloc + tlen + readlength <= rightCl:
					rsupport+=1
			line = msamlist.readline()
		tsupport=lsupport+rsupport
		ispass 		   = 0	
		if(lsupport > 0 and rsupport > 0):
			ispass =1
			passNum+=1
		else:
			reffastaleft=vcfcontent[locName][2][0:leftCl]
			if(reffastaleft.count('N')>=(len(reffastaleft)/2)):
				if(rsupport > 0):
					ispass=1
			reffastaright = vcfcontent[locName][2][rightCl:len(vcfcontent[locName][2])]
			if(reffastaright.count('N')>=(len(reffastaright)/2)):
				if(lsupport > 0):
					ispass=1
		elem = vcfcontent[locName]
		if ispass ==1:
			fil.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7};FLSUP={8};FRSUP={9};FSUP={10}\n".format(chrName, location, elem[4], elem[5], elem[6], elem[7], "PASS", elem[9], lsupport, rsupport, tsupport) )
			fil2.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6};FLSUP={7};FRSUP={8};FSUP={9}\n".format(locName, elem[4], elem[5], elem[6], elem[7], "PASS", elem[9], lsupport, rsupport, tsupport) )
			passed+=1
		else:
			fil.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7};FLSUP={8};FRSUP={9};FSUP={10}\n".format(chrName, location, elem[4], elem[5], elem[6], elem[7], "FAIL", elem[9], lsupport, rsupport, tsupport) )
			fil2.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6};FLSUP={7};FRSUP={8};FSUP={9}\n".format(locName, elem[4], elem[5], elem[6], elem[7], "FAIL", elem[9], lsupport, rsupport, tsupport) )
			failed+=1
	fil.close()
	fil2.close()
	os.system("grep PASS "+FILE+"_filtered | awk '{print $2\"\t\"$4;}' | sort -k 1,1n > "+FILE+"_filtered_PASS_loc")
	os.system("sort -k 1,1 "+FILE+"_filtered_for_setcov > "+FILE+"_filtered_for_setcov.sorted")
	#Update VCF header file with new variables obtained from filtering
	vcf_header = open("{0}/vcf_header".format(workdir),"r").readlines()
	head = open("{0}/vcf_header_f".format(workdir),"w")
	i=0
	while i < len(vcf_header):
		if vcf_header[i][2:8]=="contig":
			break
		head.write(vcf_header[i])
		i+=1
	head.write("##INFO=<ID=FLSUP,Number=1,Type=Integer,Description=\"Number of left supporting reads in filtering\">\n")
	head.write("##INFO=<ID=FLRSUP,Number=1,Type=Integer,Description=\"Number of right supporting reads in filtering\">\n")
	head.write("##INFO=<ID=FSUP,Number=1,Type=Integer,Description=\"Number of total supporting reads in filtering\">\n")
	while i < len(vcf_header):
		head.write(vcf_header[i])
		i+=1
	head.close()
#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

