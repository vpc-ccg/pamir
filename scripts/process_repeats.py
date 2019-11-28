

import sys


reps = []
prev_event = ""
content = {}
for line in sys.stdin:
    fields = line.rstrip().split("\t")
    event = fields[0]
    if event != prev_event:
        uniq_count = 0
        if len(reps) > 0:
            for prev,nxt in  zip([["-1"],*reps[:-1]],reps):
                if prev[0] == "-1":
                    uniq_count+= int(nxt[1]) - int(nxt[7])
                else:
                    uniq_count+= max(int(nxt[1]) - int(prev[2]),0)
         #       print(nxt)
            uniq_count+= int(nxt[8]) - int(nxt[2])
        #    print(event,uniq_count,sep="\t")

            prev_ins_size = int(nxt[8]) - int(nxt[7])
            print("\t".join(nxt[6:10]),end="")
            print(";unique_content={:.2f};repeat_ranges=".format(uniq_count/prev_ins_size),end="")
            types = []
            for k,v in content.items():
                types.append("{}({})".format(k,",".join(["{}-{}".format(x[0],x[1]) for x in v])))
            print(":".join(types))
        prev_event = event
        reps = [fields]
        content = {}
    else:
        reps.append(fields)
    repeat_type = fields[3].split("::")[1]
    if repeat_type not in content:
        content[repeat_type] = []
    content[repeat_type].append((int(fields[1]),int(fields[2])))
        #print(line.rstrip())

if content != {}:
    for prev,nxt in  zip([["-1"],*reps[:-1]],reps):
        if prev[0] == "-1":
            uniq_count+= int(nxt[1]) - int(nxt[7])
        else:
            uniq_count+= max(int(nxt[1]) - int(prev[2]),0)
 #       print(nxt)
    uniq_count+= int(nxt[8]) - int(nxt[2])
#    print(event,uniq_count,sep="\t")

    prev_ins_size = int(nxt[8]) - int(nxt[7])
    print("\t".join(nxt[6:10]),end="")
    print(";unique_content={:.2f};repeat_ranges=".format(uniq_count/prev_ins_size),end="")
    types = []
    for k,v in content.items():
        types.append("{}:({})".format(k,",".join(["{}-{}".format(x[0],x[1]) for x in v])))
    print(",".join(types))

