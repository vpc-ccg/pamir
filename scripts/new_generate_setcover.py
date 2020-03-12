import sys

#chr1    16739   .       G       GGGGGGGCGGCCGGGGCGGGTGGCCTCCGGGGGCGGGGGTTTGGGGGGGGGGGGGGGGTGCCCGGGCCCCCACCCAGGGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG  97.330482       PASS    Cluster=6;Support=3
def main():
    sam_file_path = sys.argv[1]
    vcf_file_path = sys.argv[2]
    out_prep_path = sys.argv[3]
    
    do = True

    events = {}

    with open(vcf_file_path, 'r') as vcf_hand:
        for line in vcf_hand:
            if line[0] == "#":
                continue
            if "PASS" not in line:
                continue

            fields = line.rstrip().split("\t")
            nom = fields[0] + "-" + fields[1]
            nom_template = nom + "-{}"
            cnt = 2
           
            if do:
                while nom in events:
                   nom = nom_template.format(cnt)
                   cnt+=1
            if nom not in events:
                events[nom] = set()

    with open(sam_file_path, 'r') as sam_hand:
        
        for line in sam_hand:
            if line[0] == "@":
                continue

            fields = line.rstrip().split("\t")
            if int(fields[8]) < 0:
                continue
            event_id = fields[2]
            if not do:
                event_id = event_id.split("-")[0] + "-" + event_id.split("-")[1]

            if event_id in events:
                events[event_id].add(fields[0])


    cluster_id = 1
    with open(out_prep_path, 'w') as hand:
        for k,v in events.items():
            if len(v) == 0:
                continue
            print(cluster_id, len(v), k, sep =" ", file=hand)
            for vv in v: #sorted(v, key = lambda x: int(x.split(":")[-1].split("/")[0])):
                print(vv, file=hand)
            cluster_id +=1

    return 0



if __name__ == "__main__":
    exit(main())
