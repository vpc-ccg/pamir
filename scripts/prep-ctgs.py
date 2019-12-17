import sys
import math
import os


sr   = open( snakemake.input[0], 'r')
sw   = open( snakemake.output[0], 'w')
sw_s = open( snakemake.output[1], 'w')
sc   = open( snakemake.output[2], 'w')

cnt_num = 1;
sw_s.write(">"+cnt_num+"\n")

line=sr.readline()
while line != '':
    cname= line.split(" ")[0]
    if cname.find('\n')==-1:
        cname=cname+'\n'

    if (cnt_loc+len(line)-1 > 1000000000):
        loc = 1;
        cnt_num ++;
        sw_s.write(">"+cnt_num+"\n")

    sc.write(cnt_num+"\t"+cname[1:len(cname)-1]+"\t"+str(cnt_loc)+'\n')

    line=sr.readline()

    sw_s.write(line)

    sw.write(cname)
    sw.write(line)

    cnt_loc+=len(line)-1
    line=sr.readline()

sr.close()
sw.close()
sw_s.close()
sc.close()

