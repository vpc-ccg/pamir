#!/usr/bin/env python
#import os, sys, errno, argparse, subprocess, fnmatch, configparser, shutil
import os, sys, errno, argparse, subprocess, fnmatch, shutil

import json

def usage():
    print('\nUsage: python filtering.py VCF REF mrsFAST-min mrsFAST-max workdir TLEN THREADS samplenum  read_len')
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
        print("Error: " +ref + " not found in reference genome")
        exit(-1)
    return ref_dict[ref][start:end]
################################
def main():
    args = sys.argv[1:]
    if len(args) != 7:
        usage()
    REF            =    sys.argv[2]
    FILE        =    sys.argv[1]
    #mrsfast min
#    MIN            =    sys.argv[3]
    #mrsfast max
#    MAX            =    sys.argv[4]
    workdir        =     sys.argv[3]
    #how many bp before the breakpoint and after the breakpoint on ref to get. (1000 for now)
    TLEN        =    int(sys.argv[4])
    THREADS        =   int(sys.argv[5])
    #samples_txt = sys.argv[8]
    all_reads_fastq = sys.argv[6]
    MRSFAST        = "mrsfast"
#    readlength = int(sys.argv[9])
    config_json = sys.argv[7]
    with open(config_json , 'r') as hand:
        cfg = json.load(hand)

    MIN = cfg["tlen_min"]
    MAX = cfg["tlen_max"]
    readlength = cfg["read_len"]

    script_folder=os.path.dirname(os.path.abspath(__file__))
    #sys.stderr.write(os.getcwd())
    RECALIBRATE = script_folder + "/../pamir recalibrate"
    folder  ="{0}/filtering".format(workdir)
    os.system("mkdir -p {0}".format(folder))
    start=1
    PREP_CTGS=script_folder + "/prep-ctgs.py"
    #f = open("{0}/allinsertions.fa".format(folder),"w")
    #coor = open("{0}/allinsertions.coor".format(folder),"w")
    #f.write(">1\n")
    prev_loc=""
    prev_chrN=""
    failed = 0
    passed = 0
    a = 2
    vcfcontent = dict()
    fil = open (FILE + "_filtered","w")
    fil2 = open (FILE + "_filtered_for_setcov","w")

    ref_dict, ref_list = load_fasta(REF)
    f_fasta = open("{0}/seq.fa".format(folder),"w")


	# todo: clean-up pack/unpack segment
	

    with open(FILE) as insertions:
        for line in insertions:
            elem_ins=line.split()
            if len(elem_ins) > 0:
                chrN    = elem_ins[0]
                loc     = elem_ins[1]
                r2        = elem_ins[2]
                r3        = elem_ins[3]
                r4         = elem_ins[4]
                r5         = elem_ins[5]
                r6         = elem_ins[6]
                r7         = elem_ins[7]
                tmplst     = r7.split(";")
#                if(len(tmplst) < 2):
#                    print(line.rstrip())
#                    print(tmplst)
#                    exit(-1)
#                length     = int(tmplst[1].split("=")[1])
                length = len(r4) - 1
#                seq        = tmplst[5].split("=")[1]
                seq = r4[1:]
                leftCl  = int(TLEN-readlength)
                rightCl = int(TLEN+length+1)
                be=int(loc)-TLEN
                if be<0:
                    be = 0
                en=int(loc)+TLEN
                left  = get_bed_seq( ref_dict, chrN, be, int(loc))
                right = get_bed_seq( ref_dict, chrN, int(loc), en)
                seqfa = left + seq + right
     #           f.write(seqfa)
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
      #          coor.write("{0}-{1}\t{2}\n".format(chrN, loc2, start ))
                f_fasta.write(">{0}-{1}\n{2}\n".format(chrN, loc2, seqfa ))
                start=start+len(left)+len(seq)+len(right)
    #f.close()
    #coor.close()
    f_fasta.close()	
    os.system("python {0} {1}/{2} {1}/{3} {4}".format( PREP_CTGS, folder, "seq.fa", "allinsertions.fa", 1000))
    
    input_file = all_reads_fastq
    

    cmd = MRSFAST + " --search {0}/allinsertions.fa --pe --min {1} --max {2} -o {0}/seq.mrsfast.sam -e 3 --threads {4} --disable-sam-header --disable-nohits --seq {3} > {0}/seq.mrsfast.sam.log".format(folder, MIN, MAX,input_file,THREADS)
    #print(cmd,file=sys.stderr)
    sys.stderr.write(os.getcwd()+"\n")
    sys.stderr.write(script_folder + __file__ +"\n" + "\n")
 
    sys.stderr.write( cmd + "\n" + "\n")
    os.system(MRSFAST + " --index {0}/allinsertions.fa > {0}/mrsfast.index.log".format(folder))

#    with open(samples_txt,'r') as hand:
#        for line in hand:
#            input_file = input_file + "{0}/{1}/{1}.all_interleaved.fastq".format(workdir, line.rstrip())





    fname = "{0}/run_filtering.sh".format(workdir)
    f     = open(fname,"w")
    f.write("#!/bin/bash\n")
    f.write(cmd+'\n')
    f.close()
    cmd           =  'cat {0} | xargs -I CMD --max-procs=1 bash -c CMD'.format(fname)    
    os.system(cmd)
    os.system("{1} {0}/allinsertions.fa.lookup {0}/seq.mrsfast.sam {0}/seq.mrsfast.recal.sam".format(folder, RECALIBRATE) )
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
        contig_parts = locName.split("-")

        chrName = contig_parts[0]
        location = contig_parts[1]
        if contig_parts[0] == "HLA":
            location = contig_parts[2]
            chrName = "{}-{}".format(contig_parts[0],contig_parts[1])
#        first_sep = locName.find("-")
#        last_sep = locName.rfind("-")
#        chrName = locName[0:first_sep]
#        if(first_sep ==last_sep):
#            last_sep = len(locName)
#        location = locName[first_sep+1:last_sep]
        firstloc = int(splitmsam[3])
        tlen     = int(splitmsam[8])
        rightCl = int(TLEN+int(vcfcontent[locName][0])+1)
        if flag & 2 == 2:
            if firstloc < TLEN - 10 and firstloc + readlength >= TLEN + 10:
                lsupport+=1
            if firstloc < rightCl - 10 and firstloc + readlength > rightCl + 10:
                rsupport+=1
        line = msamlist.readline()
        while(line !=''):
            splitmsam    = line.split()
            flag         = int(splitmsam[1])
            nextlocName = splitmsam[2]
            nextloc = int(splitmsam[3])
            tlen     = int(splitmsam[8])
            if(nextlocName != locName):
                break;
            if flag & 2 ==2:
                #if nextloc <= leftCl and nextloc + tlen >= (TLEN+1):
                if nextloc < TLEN - 10 and nextloc + readlength >= TLEN + 10:
                    lsupport+=1
                #if nextloc >= rightCl and nextloc + tlen + readlength <= rightCl:
                if nextloc < rightCl - 10 and nextloc + readlength > rightCl + 10:
                    rsupport+=1
            line = msamlist.readline()
        tsupport=lsupport+rsupport
        ispass            = 0    
        if(lsupport > 0 and rsupport > 0):
            ispass =1
            passNum+=1

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
    return 0
#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

