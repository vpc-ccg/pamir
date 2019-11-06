#!/usr/bin/env python
#-*- encoding:utf-8 -*-

import sys, os, math, operator

##########
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


##########
def usage():
    print('Usage: python '  + sys.argv[0] + ' input.vcf output.vcf rmdup')
    print('\t' + sys.argv[0] + ': sorting VCF by breakpoints, insertion lengths and insertion contents')
    print('\tinput: vcf file by pamir')
    print('\toutput: sorted by chromosome ')
    print('\trmdup: set 1 if you also want to remove the duplicated insertions (same loc, length, sequence), otherwise 0.')
    sys.exit(-1)

##########
def parse_ins_vcf( info_str ):
    sv_length = 0 
    sv_seq = ""
    info_list = info_str.split(";")
    for item in info_list:
        tmp_list = item.split("=")
        if ( "SVLEN" == tmp_list[0]):
            sv_length = int(tmp_list[1])
        elif ("SEQ" == tmp_list[0]):
            sv_seq = tmp_list[1]

    return sv_length, sv_seq

##########
def main():
    args = sys.argv[1:]
    if len(args) != 3:
        usage()

    
    chr_list = []
    cur_chr= ""
    vcf_dict = {}
    sw = open( sys.argv[2], 'w')
    sr = open( sys.argv[1], 'r')
    rmdup = int(sys.argv[3])
    for line in sr:
        if ("#" == line[0]):
            sw.write(line)
            continue
        
        t_list = line.strip().split("\t")
        chr_str = t_list[0]
        # Get chrmosome order in original vcf file 
        if ( chr_str != cur_chr):
            chr_list.append( chr_str)
            cur_chr = chr_str

        bp = int(t_list[1])
#        sv_length, sv_seq = parse_ins_vcf( t_list[7] )
        sv_seq = t_list[4][1:]
        sv_length = len(sv_seq)
        sv_token=[ bp, sv_length, sv_seq, line ]

        if chr_str not in vcf_dict:
            vcf_dict[ chr_str ] = []

        vcf_dict[ chr_str ].append( sv_token )


    #for x,y in vcf_dict.iteritems():
    for x in chr_list:
        y = vcf_dict[x]
        y.sort(key = operator.itemgetter(0, 1, 2) )
        if rmdup == 0 :
            for item in y:
                sw.write( item[3] )
        elif rmdup == 1:
            i = 0
            while i+1 < len(y):
                if not(y[i][0] ==y[i+1][0] and y[i][1]==y[i+1][1] and y[i][2]==y[i+1][2]):
                    sw.write(y[i][3])
                i+=1
            if not(y[i-1][0] ==y[i][0] and y[i-1][1]==y[i][1] and y[i-1][2]==y[i][2]):
                sw.write(y[i][3])


    sr.close()
    sw.close()

if __name__ == "__main__":
    sys.exit(main())
