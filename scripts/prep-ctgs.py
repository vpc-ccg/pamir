import sys, math, os

#input fasta spacer_size

origin_fasta = open(sys.argv[1] , 'r')
merged_fasta = open(sys.argv[2], 'w')
merged_fasta_idx = open(sys.argv[2]+".lookup", 'w')
spacer_size=int(sys.argv[3])

cnt_num = 1;
cnt_loc=1;
content_rem = ""

merged_fasta.write(">"+str(cnt_num)+"\n")

# reading the first line of the fasta file
line=origin_fasta.readline()
spacer=spacer_size * "N"

while line != '':

    #extracting the name up to first white space; ignoring the comments
    cname= line.split(" ")[0]
    if cname.find('\n')==-1:
        cname=cname+'\n'


    # reading the input until start of the next contig or end of file
    content = ""
    line=origin_fasta.readline()
    while (line != "" and line.find(">") == -1):
        content = content + line.replace("\n","")
        line=origin_fasta.readline()
    content = content + spacer

    # starting a new chromosome
    if (cnt_loc-1 +len(content) > 1000000000):
        merged_fasta.write(content_rem+"\n")
        cnt_num = cnt_num+1;
        merged_fasta.write(">"+str(cnt_num)+"\n")
        cnt_loc = 1
        content_rem = ""


    #formatted dump
    dump_content = content_rem + content
    dump_len = len (dump_content)
    if ( dump_len > 80 ):
        cut = (int)(dump_len / 80)*80;
        content_rem = dump_content [cut:]
        [ merged_fasta.write( dump_content[i:i+80] + "\n") for i in range(0, cut, 80)]
    else:
        content_rem = dump_content

    merged_fasta_idx.write("{}\t{}\t{}\n".format(cnt_num, cname[1:len(cname)-1], cnt_loc))

    cnt_loc+=len(content)

# Dump the final content out
if (content_rem != ""):
    merged_fasta.write(content_rem+"\n")

origin_fasta.close()
merged_fasta.close()
merged_fasta_idx.close()

