#!/usr/bin/env python
import os, sys, errno, argparse, subprocess, fnmatch, ConfigParser, shutil
def main():
	REF=sys.argv[2]
	FILE=sys.argv[1]
	if(os.path.isfile(FILE+"_filtered")):
		os.unlink(FILE+"_filtered")
	with open(FILE) as insertions:
#		count=0
		fil = open (FILE + "_filtered","a")
		passed=0
		failed=0
		for line in insertions:
			elem_ins=line.split()
			if len(elem_ins) > 0:
				chrN    = elem_ins[0]
				loc     = elem_ins[1]
				length  = elem_ins[2]
				seq     = elem_ins[3]
				folder  ="fpcheck/{0}".format(loc)
				leftCl  = int(900)
				rightCl = int(1000+int(length)+1)
				os.system("mkdir -p {0}".format(folder))
				os.system("echo {0} > {1}/seq.txt".format(seq, folder))
				be=int(loc)-1-1001
				en=int(loc)+1000
				open("{0}/left.bed".format(folder),"w").write("{0}\t{1}\t{2}\n".format(chrN,be,int(loc)-1))
				open("{0}/right.bed".format(folder),"w").write("{0}\t{1}\t{2}\n".format(chrN,int(loc)-1,en))
				os.system("samtools faidx {0}".format(REF))
				os.system("bedtools getfasta -bed {0}/left.bed -fi {1} -fo {0}/left.fa".format(folder, REF))
				os.system("bedtools getfasta -bed {0}/right.bed -fi {1} -fo {0}/right.fa".format(folder, REF))
				left=open("{0}/left.fa".format(folder),"r").readlines()[-1]
				right=open("{0}/right.fa".format(folder),"r").readlines()[-1]
				open("{0}/newgen.fa".format(folder),"w").write(">genome_{0}\n{1}{2}{3}".format(loc, left[:len(left)-1], seq, right[:len(right)-1] ))
				os.system("../mrsfast --index {0}/newgen.fa > {0}/mrsfast.index.log".format(folder))
				os.system("samtools faidx {0}/newgen.fa".format(folder))
#				os.system("bwa index {0}/newgen.fa".format(folder))	
				os.system("../mrsfast --search {0}/newgen.fa --pe --min 300 --max 500 -n 1 -o {0}/seq.mrsfast.sam -e 3 --seq all_interleaved.fastq --disable-nohits > {0}/.seq.mrsfast.sam.log".format(folder))
#				os.system("bwa mem {0}/newgen.fa all_1.fastq all_2.fastq > {0}/seq.bwa.sam 2> {0}/seq.bwa.sam.log".format(folder))
				msamlist = open("{0}/seq.mrsfast.sam".format(folder),"r").readlines()
#				msamlist = open("{0}/seq.bwa.sam".format(folder),"r").readlines()
				lsupport=0
				rsupport=0
				tsupport=0
				if(len(msamlist)>2):
					os.system("../sniper sort {0}/seq.mrsfast.sam {0}/seq.mrsfast.sam.sorted".format(folder))
#					fbwa= open("{0}/seq.bwa.sam_wonohit".format(folder),"w") 
#					for elem in range(2,len(msamlist)):
#						if(msamlist[elem].split()[2]!="*"):
#							fbwa.write(msamlist[elem])
#					fbwa.close()
#					os.system("../sniper sort {0}/seq.bwa.sam_wonohit {0}/seq.bwa.sam.sorted".format(folder))

					msamlist = open("{0}/seq.mrsfast.sam.sorted".format(folder),"r").readlines()
#					msamlist = open("{0}/seq.bwa.sam.sorted".format(folder),"r").readlines()
					firstloc       = msamlist[0].split()[3]
					lastloc 	   = msamlist[len(msamlist)-1].split()[3]
					ispass 		   = 0	
					if(int(firstloc)<=leftCl and int(lastloc)>=rightCl):
						ispass =1
					else:
						reffastaleft=""
						reflinesleft = open("{0}/left.fa".format(folder),"r").readlines()
						for p in range (1,len(reflinesleft)):
							reffastaleft += reflinesleft[p][:len(reflinesleft[p])-1]
						if(reffastaleft.count('N')>=(len(reffastaleft)/2)):
							if(int(firstloc)<=leftCl or int(lastloc)>=rightCl):
								ispass=1
						reffastaright =""
						reflinesright = open("{0}/right.fa".format(folder),"r").readlines()
						for p in range (1,len(reflinesright)):
							reffastaright += reflinesright[p][:len(reflinesright[p])-1]
						if(reffastaright.count('N')>=(len(reffastaright)/2)):
							if(int(firstloc)<=leftCl or int(lastloc)>=rightCl):
								ispass=1
					if ispass ==1:
						for a in range (0,len(msamlist)):
							l = int(msamlist[a].split()[3])
							if l <=leftCl:
								lsupport+=1
							if l>=rightCl:
								rsupport+=1
						tsupport=lsupport+rsupport
						fil.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tPASS\n".format(chrN, loc, seq, length, lsupport, rsupport, tsupport) )
						passed+=1
					else:
						fil.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tFAIL\n".format(chrN, loc, seq, length, lsupport, rsupport, tsupport) )
						failed+=1
						shutil.rmtree(folder,ignore_errors=True)
				else:
					fil.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tFAIL\n".format(chrN, loc, seq, length, lsupport, rsupport, tsupport) )
					failed+=1
					shutil.rmtree(folder,ignore_errors=True)
		fil.close()
		os.system("grep PASS sniper_part_updated.vcf_sorted_wodups_filtered | awk '{print $2\"\t\"$4;}' | sort -k 1,1n > sniper_part_updated.vcf_sorted_wodups_filtered_PASS_loc")

#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

