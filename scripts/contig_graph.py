import sys

CIGAR_STR = "DHIMNPSX="
def check_cchar(ch):
    for c in CIGAR_STR:
        if ch == c:
            return ch

    return 0
def find_end(pos,cgr):

    prv = -1
    pos = 0
    end = pos
    while pos < len(cgr):
        ch = check_cchar(cgr[pos])

        if ch != 0:
            ln = int(cgr[1+prv:pos])
            if ch == 'M' or ch == 'D' or ch == '=' or ch =='X':
                end += ln
            prv = pos
        pos+=1
    return end
reads = {}
def add_read(reads, line):
    fields = line.split("\t")
    rid = fields[0]
    if rid not in reads:
        reads[rid] = ([],[])
    flag = int(fields[1])


    is_first  =  (flag & 64) != 0
    is_second =  (flag & 128) != 0
    if is_first:
        reads[rid][0].append(fields)
    if is_second:
        reads[rid][1].append(fields)

for line in sys.stdin:
    add_read(reads,line)

rg = {}

def add_edge(f,s):
    if f[2] not in rg:
        rg[f[2]] ={}
    if s[2] not in rg:
        rg[s[2]] ={}
    

    f_end = find_end(int(f[3]),f[5])
    s_end = find_end(int(s[3]),s[5])
    if s[2] not in rg[f[2]]:
        rg[f[2]][s[2]] = (int(f[3]),int(f[3])+f_end,1)
    else:
        prev=rg[f[2]][s[2]]
        rg[f[2]][s[2]]=(min(prev[0],int(f[3])),max(prev[1],int(f[3])+f_end),prev[2]+1)

    if f[2] not in rg[s[2]]:
        rg[s[2]][f[2]] = (int(s[3]),int(s[3])+s_end,1)
    else:
        prev=rg[s[2]][f[2]]
        rg[s[2]][f[2]]=(min(prev[0],int(s[3])),max(prev[1],int(s[3])+s_end),prev[2]+1)

for r,a in reads.items():
    for f in a[0]:
        for s in a[1]:
            if f[2] != s[2]:
                add_edge(f,s)

def visit(graph, visited, node, target, depth):
    if node == target:
        return depth
    for v in graph[node]:
        if v not in visited:
            ret = visit(graph,visited,v,target,depth+1)
            if ret > 0:
                return ret
    return -1
def bfs(graph,n1,n2):
    visited = {n1:1}
    for v in graph[n1]:
        if v not in visited:
            ret = visit(graph,visited,v,n2,1)
            if ret > 0:
                return ret
    return -1
#for k,v in rg.items():
#    for a,b in v.items():
#        print("{}\t{}\t{}\t{}-{}".format(k,a,b[2],b[0],b[1]))


with open(sys.argv[1],'r') as hand:
    for line in hand:
        fields = line.rstrip().split("\t")
        dist = bfs(rg,fields[0],fields[1])

        print("{}\t{}".format(line.rstrip(),dist))
