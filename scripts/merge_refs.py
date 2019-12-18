

import sys

#Does not like multiline fastas
lookup = {}
line = sys.stdin.readline()
index_path = sys.argv[1]
while line:
    fpath = line.rstrip()
    with open(fpath, 'r') as hand:
        fline = hand.readline()
        while(fline):
            cid = fline.rstrip()
            seq = hand.readline()
            if cid not in lookup:
                lookup[cid] = []
                print(cid)
                print(seq,end="")
            lookup[cid].append(fpath[fpath.rfind("/"):fpath.rfind(".")])
            fline = hand.readline()
    line = sys.stdin.readline()
with open(index_path, 'w') as whand:
    for k,v in lookup.items():
        print(k,end="",file=whand)
        for vv in v:
            print("\t{}".format(vv),end="",file=whand)
        print("",file=whand)
