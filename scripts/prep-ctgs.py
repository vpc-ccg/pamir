import sys
import math
import os



#sr   = open( "deneme.fa", 'r')
#sw   = open( "output1", 'w')
#sw_s = open( "output2", 'w')
#sc   = open( "output3", 'w')

sr   = open( snakemake.input[0], 'r')
sw   = open( snakemake.output[0], 'w')
sw_s = open( snakemake.output[1], 'w')
sc   = open( snakemake.output[2], 'w')

cnt_num = 1;
cnt_loc=1;
sw_s.write(">"+str(cnt_num)+"\n")

line=sr.readline()
while line != '':
    cname= line.split(" ")[0]
    if cname.find('\n')==-1:
        cname=cname+'\n'

    line=sr.readline()
    if (cnt_loc+len(line)-1 > 1000000000):
        cnt_loc = 1;
        cnt_num = cnt_num+1;
        sw_s.write(">"+str(cnt_num)+"\n")

    sw_s.write(line+"NNNNNNNNNN\n")
    sc.write(str(cnt_num)+"\t"+cname[1:len(cname)-1]+"\t"+str(cnt_loc)+'\n')



    sw.write(cname)
    sw.write(line)

    cnt_loc+=len(line)-1 + 10 ## for Ns
    line=sr.readline()

sr.close()
sw.close()
sw_s.close()
sc.close()

