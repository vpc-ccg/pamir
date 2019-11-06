#!/usr/bin/env python
import os, sys, errno, argparse, subprocess, fnmatch, configparser, shutil

def usage():
    print('\nUsage: python genotyping.py VCF REF SEQ1.fastq SEQ2.fastq SAMPLENAME mrsFAST-min mrsFAST-max workdir TLEN THREADS read_len mrsfast recalibrate ')
    sys.exit(-1)
################################
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
# Get sequences from ref in a bed-fashion
#########################
def get_bed_seq( ref_dict, ref, start, end):
    if ref not in ref_dict:
        print("Error: " +ref + " not found in reference genome")
        exit(-1)
    return ref_dict[ref][start:end]
#########################
def main():
    args = sys.argv[1:]
    if len(args) !=13:
        usage()
    REF=sys.argv[2]
    #FILE needs to be sorted according to 1st and 2nd columns
    FILE=sys.argv[1]
    SEQ1=sys.argv[3]
    SEQ2=sys.argv[4]

    EXT=sys.argv[5]
    MIN=sys.argv[6]
    MAX=sys.argv[7]
    workdir=sys.argv[8]
    TLEN = int(sys.argv[9])
    THREADS = int(sys.argv[10])
    readlen = int(sys.argv[11])
    MRSFAST= sys.argv[12]
    RECALIBRATE = sys.argv[13]
    folder  ="{0}/genotype".format(workdir)
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
    ref_dict, ref_list = load_fasta(REF)

    fil = open ("{0}/insertions.out_wodups_filtered_setcov".format(workdir) + "_genotype_"+EXT,"w")
    with open(FILE) as insertions:
        for line in insertions:
            elem_ins=line.split()
            if len(elem_ins) > 0:
                chrN        = elem_ins[0]
                loc         = elem_ins[1]
                tmplist        = elem_ins[7].split(";")

                #seq         = tmplist[5].split("=")[1]
                seq = elem_ins[4][1:]
                length = len(seq)
                #length        = tmplist[1].split("=")[1]
                

                be=int(loc)-TLEN
                if be<0:
                    be =0
                en=int(loc)+TLEN
                left  = get_bed_seq( ref_dict, chrN, be, int(loc))
                right = get_bed_seq( ref_dict, chrN, int(loc), en)
                seqfaref = left + right
                seqfains = left + seq + right
                f.write(seqfains)
                f2.write(seqfaref)
                loc2 =loc
                if(chrN==prev_chrN and loc==prev_loc):
                    loc2 = loc+"-"+str(a)
                    key = chrN + "-" + loc2
                    vcfcontent[key]=[]
                    vcfcontent[key].append(elem_ins[2])
                    vcfcontent[key].append(elem_ins[3])
                    vcfcontent[key].append(elem_ins[4])
                    vcfcontent[key].append(elem_ins[5])
                    vcfcontent[key].append(elem_ins[6])
                    vcfcontent[key].append(elem_ins[7])
                    vcfcontent[key].append(0)
                    vcfcontent[key].append(0)
                    a+=1
                else:
                    prev_loc =loc
                    prev_chrN =chrN
                    key = chrN + "-" + loc
                    vcfcontent[key]=[]
                    vcfcontent[key].append(elem_ins[2])
                    vcfcontent[key].append(elem_ins[3])
                    vcfcontent[key].append(elem_ins[4])
                    vcfcontent[key].append(elem_ins[5])
                    vcfcontent[key].append(elem_ins[6])
                    vcfcontent[key].append(elem_ins[7])
                    vcfcontent[key].append(0)
                    vcfcontent[key].append(0)
                    a=2
                coor.write("{0}-{1}\t{2}\n".format(chrN, loc2, start ))
                coor2.write("{0}-{1}\t{2}\n".format(chrN, loc2, start2 ))
                start=start+len(left)+len(seq)+len(right)
                start2=start2+len(left)+len(right)
    f.close()
    f2.close()
    coor.close()
    coor2.close()
    os.system(MRSFAST + " --index {0}/allref.fa > {0}/mrsfast.index.log".format(folder))
    os.system(MRSFAST + " --search {0}/allref.fa --seqcomp --pe --min {4} --max {5} -n 50 --threads {6} -o {0}/seq.mrsfast.ref.{1}.sam -e 3 --seq1 {2} --seq2 {3} --disable-sam-header --disable-nohits > {0}/seq.mrsfast.ref.{1}.sam.log".format(folder, EXT, SEQ1, SEQ2, MIN, MAX, THREADS))
    os.system(RECALIBRATE + " {0}/allref.coor {0}/seq.mrsfast.ref.{1}.sam {0}/seq.mrsfast.ref.{1}.recal.sam".format(folder,EXT))
    os.system("sort -S 500G -T /cs/compbio2/ -k 3,3 -k 4,4n {0}/seq.mrsfast.ref.{1}.recal.sam > {0}/seq.mrsfast.ref.{1}.recal.sam.sorted".format(folder,EXT))
    msamlist = open("{0}/seq.mrsfast.ref.{1}.recal.sam.sorted".format(folder,EXT),"r")

    os.system(MRSFAST + " --index {0}/allinsertions.fa > {0}/mrsfast.index2.log".format(folder))
    os.system(MRSFAST + " --search {0}/allinsertions.fa --seqcomp --pe --min {4} --max {5} -n 50 --threads {6} -o {0}/seq.mrsfast.ins.{1}.sam -e 3 --seq1 {2} --seq2 {3} --disable-sam-header --disable-nohits > {0}/seq.mrsfast.ins.{1}.sam.log".format(folder, EXT, SEQ1, SEQ2, MIN, MAX, THREADS))
    os.system(RECALIBRATE + " {0}/allinsertions.coor {0}/seq.mrsfast.ins.{1}.sam {0}/seq.mrsfast.ins.{1}.recal.sam".format(folder,EXT))
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
        first_sep = locName.find("-")
        last_sep = locName.rfind("-")
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
        vcfcontent[locName][6]=refsupport
    line2 = msamlist2.readline()

    while (line2 !=''):
        altsupport_left=0
        altsupport_right=0
        splitline2 = line2.split()
        flag2 = int(splitline[1])
        locName2 = splitline2[2]
        tmp = vcfcontent[locName2][5].split(";")
        tmplist = tmp[1].split("=")
        lenins = tmplist[1]
        secondbreakpoint = breakpoint+int(lenins)
        first_sep2 = locName2.find("-")
        last_sep2 = locName2.rfind("-")
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
        vcfcontent[locName2][7]= altsupport
    for a in vcfcontent:
        elem = vcfcontent[a]
        chrNameLoc=a.split('-')
        if chrNameLoc[0] == "HLA":
            this_chr = "-".join([chrNameLoc[0],chrNameLoc[1]])
            this_pos = chrNameLoc[2]
        else:
            this_chr = chrNameLoc[0]
            this_pos = chrNameLoc[1]
        ratio=0
        if elem[6]==0 and elem[7]==0:
            gtype="0/0"
        else:
            ratio = (float(elem[7])-float(elem[6]))/(float(elem[7])+float(elem[6]))
            if(ratio >= 0.3):
                gtype="1/1"
            elif(ratio <=-0.3):
                gtype="0/0"
            else:
                gtype="0/1"
        fil.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7};SAMPLE={12};GRSUP={8};GISUP={9};GRATIO={10}\tGT\t{11}\n".format(this_chr, this_pos, elem[0], elem[1], elem[2], elem[3], elem[4], elem[5], elem[6], elem[7],str(ratio),gtype,EXT) )
    fil.close()

    #Update VCF header file with new variables obtained from filtering

    os.system("cat {0}/insertions.out_wodups_filtered_setcov_genotype_{1} > {0}/insertions_genotype_{1}.vcf.nohead".format(workdir, EXT))
#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

