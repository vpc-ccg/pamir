import sys
import math
import os


sr   = open( snakemake.input[0], 'r')
sw   = open( snakemake.output[0], 'w')
sw_s = open( snakemake.output[1], 'w')
sc   = open( snakemake.output[2], 'w')
sw_s.write(">1\n")
cloc = 1

line=sr.readline()
while line != '':
    cname= line.split(" ")[0]
    if cname.find('\n')==-1:
        cname=cname+'\n'
    sc.write(cname[1:len(cname)-1]+"\t"+str(cloc)+'\n')
    sw.write(cname)
    line=sr.readline()
    cloc+=len(line)-1
    sw.write(line)
    sw_s.write(line)
    line=sr.readline()

sr.close()
sw.close()
sw_s.close()
sc.close()

