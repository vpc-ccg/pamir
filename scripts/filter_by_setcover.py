#!/usr/bin/env python
#-*- encoding:utf-8 -*-
import sys
import math
import os
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg
def usage():
    print('\nUsage: python filter_by_setcover.py fromSETCOVER VCF VCFafterSETCOVER')
    sys.exit(-1)

def main():
    args = sys.argv[1:]
    if len(args) !=3:
        usage()

    fromSetCover = open(sys.argv[1],"r")
    vcffile = open(sys.argv[2],"r")
    fout = open(sys.argv[3],"w")
    dictvcf = dict()
    for line in vcffile:
        if line[0] == "#":
            fout.write(line)
            continue
        spline = line.split()

        key = spline[0]+"-"+spline[1]
        if key in dictvcf:
            dictvcf[key].append(spline[2])
            dictvcf[key].append(spline[3])
            dictvcf[key].append(spline[4])
            dictvcf[key].append(spline[5])
            dictvcf[key].append(spline[6])
            dictvcf[key].append(spline[7])
        else:        
            dictvcf[key]=[]
            dictvcf[key].append(spline[2])
            dictvcf[key].append(spline[3])
            dictvcf[key].append(spline[4])
            dictvcf[key].append(spline[5])
            dictvcf[key].append(spline[6])
            dictvcf[key].append(spline[7])
    for line in fromSetCover:
        spline = line.split()
        if spline[0]!="Removed:":
            nameloc = spline[2].split("-")
            if(nameloc[0] == "HLA"):
                chrname = "-".join([nameloc[0],nameloc[1]])
                chrloc = nameloc[2]
            else:
                chrname = nameloc[0]
                chrloc  = nameloc[1]
            lenkeyelem = 0
            roundnum =0
            while lenkeyelem < len(dictvcf[spline[2]])-1:
                fout.write(chrname+"\t"+chrloc)
                for i in range(0,6):
                    fout.write("\t"+dictvcf[spline[2]][(roundnum*6)+i])
                    lenkeyelem+=1
                fout.write("\n")
                roundnum+=1
    fout.close()

if __name__ == "__main__":
    sys.exit(main())
