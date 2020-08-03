
def cfg_default( key, val):
    if key not in config:
        config[key] = val

def cfg_mandatory( key):   
    assert key in config, "{} is a mandatory field in config".format(key)
    
def cfg_optional( key):
    pass

def cfg_internal( key):
    assert key not in config, "{} is an internal parameter! Do not use.".format(key)

if "DEBUG" in config:
    def temp(st):
        return st

new_prep = True
assembler_binaries = {"minia":"minia","abyss":"abyss-pe","spades":"spades.py", "present":"ls"}

cfg_mandatory("path")
cfg_mandatory("reference")
cfg_mandatory("raw-data")


cfg_mandatory("input")
cfg_mandatory("population")



cfg_default("pamir_partition_per_thread",1000)
cfg_default("analysis-base","/analysis")
cfg_default("results-base","/results")

cfg_default("analysis","{}{}/{}".format(config["path"],config["analysis-base"],config["population"]))
cfg_default("results","{}{}/{}".format(config["path"],config["results-base"],config["population"]))

cfg_default("genotyping_flank_size",1000)
cfg_default("assembler","minia")
cfg_default("assembler_k", 64)
cfg_default("assembly_threads",62)
cfg_default("aligner_threads", 16)
cfg_default("other_threads",16)
cfg_default("repeatmasker_threads",63)
cfg_default("minia_min_abundance",5)
cfg_default("n_std_dev",3)
cfg_default("seq_ext", "fq.gz")

cfg_optional("pamir_min_contig_len")
cfg_optional("blastdb")
cfg_optional("DEBUG")
cfg_optional("centromeres")

sext = config["seq_ext"]
mrsfast_comp_param = "" if ".gz" not in sext else " --seqcomp "

path_names = {
"linkeddata"  : "linked-data",
"remcor"      : "pamir-remove-concordants",
"extass"      :  config["assembler"],
"blast"       : "blast-filtering",
"lenfilter"   : "length-filtering",
"pamoea"      : "pamir-oea-processing",
"partition"   : "pamir-partition",
"pamass"      : "pamir-assembly",
"pamann"      : "pamir-filtering",
"pamscv"      : "pamir-setcover",
"pamref"      : "pamir-flank-reference",
"pamgen"      : "pamir-genotyping",
}
if "SCPA" in config:
    SCPA = config["SCPA"]
else:
    SCPA = "."
index = 0
for k,v in path_names.items():
    path_names[k] = config["analysis"] + "/{:0>3d}-{}".format(index,v)
    index+=1

def tool_exists(name):
    from shutil import which
    return which(name) is not None

def assert_tool(tool):
    assert tool_exists(tool), tool + " not found in the PATH"

assert_tool(assembler_binaries[config["assembler"]])

def get_cram_name(wildcards):
    orig_cram_name = config["input"][wildcards.sample][0]
    extension = os.path.splitext(orig_cram_name)[1]
   
    if extension not in [".bam", ".cram"]:
        print("Input should be in bam or cram format.", file=sys.stderr)
        sys.exit(1)
    
    return path_names["linkeddata"] +"/"+ wildcards.sample + extension


def get_cram_index(wildcards):
    cram_name = get_cram_name(wildcards)
    extension = os.path.splitext(cram_name)[1]
    extension = {".bam":".bai",".cram":".crai"}[extension]

    return "{}{}".format(cram_name,extension)



rule make_results:
    input:
        html=config["results"]+"/index.html",
        data_js=config["results"]+"/data.js",
        summary_js=config["results"]+"/summary.js",
        bams=expand("{results}/ind/{sample}/events.bam",results=config["results"],sample=config["input"].keys()),
        bais=expand("{results}/ind/{sample}/events.bam.bai",results=config["results"],sample=config["input"].keys()),
        bed=expand("{results}/ind/{sample}/events.bed",results=config["results"],sample=config["input"].keys()),
        vcf=expand("{results}/ind/{sample}/events.vcf",results=config["results"],sample=config["input"].keys()),
        fa=expand("{results}/events.fa",results=config["results"]),
        fai=expand("{results}/events.fa.fai",results=config["results"]),
        rm=expand("{results}/events.repeat.bed",results=config["results"]),
        seq_stats=config["results"]+"/stats.js"
    
rule print_event_xml:
    input:
        bam=path_names["pamgen"]+"/{sample}/final.bam",
        fasta=path_names["pamref"]+"/events_ref.fa",
        bed  =path_names["pamgen"]+"/{sample}/events_ins_pos.bed",
    output:
        dne =path_names["pamgen"]+"/{sample}/cfg.xml"
    run:
        with open(output.dne,'w') as out_hand:
            genome_str = '<Global genome="{}" version="3">'
            bam_str='<Resource name="Reads" path="{}"/>'
            bed_str='<Resource name="Insertion Coordinates" path="{}"/>'

            covtrack_str='<Track autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" color="175,175,175" colorScale="ContinuousColorScale;0.0;19.0;255,255,255;175,175,175" fontSize="10" id="{}_coverage" name="Reads Coverage" snpThreshold="0.2" visible="true">'
            seqtrack_str='<Track clazz="org.broad.igv.sam.AlignmentTrack" id="{}" name="Reads" visible="true" displayMode="SQUISHED" >'
            bedtrack_str='<Track clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;1.0;255,255,255;0,0,178" fontSize="10" id="{}" name="Insertion Coordinates" visible="true"/>'
            xml_lines = [
                '<?xml version="1.0" encoding="UTF-8"?>',
                genome_str.format(input.fasta),
                '<Resources>',
                bam_str.format(input.bam),
                bed_str.format(input.bed),
                '</Resources>',
                '<Panel name="SeqPanel">',
                covtrack_str.format(input.bam),
                '<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="20.0" minimum="0.0" type="LINEAR"/>',
                '</Track>',
                seqtrack_str.format(input.bam),
                '<RenderOptions viewPairs="true"/>',
                '</Track>',
                bedtrack_str.format(input.bed),
                '</Panel>',
               # '<Panel name="RefPanel">',
               # '<Track clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" visible="true"/>',
               # bedtrack_str.format(input.bed),
               # '</Panel>',
                '<PanelLayout dividerFractions="0.006984866123399301,0.7558323632130385"/>',
                '</Global>'
            ]
            tabc=0;
            tab="\t";

            for line in xml_lines:

                if "</" in line:
                    tabc-=1
                print(tab*tabc,end="",file=out_hand)
                print(line,file=out_hand)
                if "<?" in line:
                    pass
                else:
                    if "<" in line and "/>" not in line and "</" not in line:
                        tabc+=1
rule make_html:
    output:
        config["results"]+"/index.html",
    params:
        scpa=SCPA,
    shell:
        "cp {params.scpa}/scripts/index.html {output}"

rule move_bams:
    input:
        path_names["pamgen"]+"/{sample}/final.bam" 
    output:
        config["results"]+"/ind/{sample}/events.bam" 
    shell:
        "cp {input} {output}"

rule move_vcf:
    input:
        rem=config["results"] + "/rem.done",
        vcf=path_names["pamgen"]+"/{sample}/final.vcf" 
    output:
        config["results"]+"/ind/{sample}/events.vcf" 
    shell:
        "cp {input.vcf} {output}"


rule move_bed:
    input:
        path_names["pamref"]+"/ind/{sample}/events_ins_pos.bed" 
    output:
        config["results"]+"/ind/{sample}/events.bed" 
    shell:
        "cp {input} {output}"

rule move_fa:
    input:
        path_names["pamref"]+"/events_ref.fa" 
    output:
        config["results"]+"/events.fa" 
    shell:
        "cp {input} {output}"

rule move_fai:
    input:
        path_names["pamref"]+"/events_ref.fa.fai" 
    output:
        config["results"]+"/events.fa.fai" 
    shell:
        "cp {input} {output}"

rule move_repeat_bed:
    input:
        path_names["pamref"]+"/events_ref.repeat.bed",
    output:
        config["results"]+"/events.repeat.bed" 
    shell:
        "cp {input} {output}"

#chr13-69551766-747f     1000    1050    Cluster=1087545;Support=4;FLSUP=6;FRSUP=8;FSUP=14       chr13-69551766-747f     1006    1084    (TA)n::Simple_repeat    0       +

rule merge_repeat_beds:
    input:
        uniq=path_names["pamref"]+"/ind/{sample}/uniq_seq_ins.bed",
        rep=path_names["pamref"]+"/ind/{sample}/repeat_overlap_ins.bed",
    output:
        temp(path_names["pamref"]+"/ind/{sample}/events_ins_pos_and_repeat.bed"),
    shell:
        "cat {input.uniq} {input.rep} | bedtools sort  -i stdin > {output}"

rule non_repeat_bed:
    input:
        bed=path_names["pamref"]+"/ind/{sample}/events_ins_pos.bed",
        rpt=path_names["pamref"]+"/events_ref.repeat.bed",
    output:
        temp(path_names["pamref"]+"/ind/{sample}/uniq_seq_ins.bed"),
    params:
        scpa=SCPA,
    shell:
        "bedtools intersect -b {input.rpt} -a {input.bed} -v  | python {params.scpa}/scripts/process_unique.py > {output}"

rule repeat_bed:
    input:
        bed=path_names["pamref"]+"/ind/{sample}/events_ins_pos.bed",
        rpt=path_names["pamref"]+"/events_ref.repeat.bed",
    output:
        temp(path_names["pamref"]+"/ind/{sample}/repeat_overlap_ins.bed"),
    params:
        scpa=SCPA,
    shell:
        "bedtools intersect -a {input.rpt} -b {input.bed} -wb | python {params.scpa}/scripts/process_repeats.py > {output}"


rule make_unique_repeats:
    input:
        path_names["pamref"]+"/events_ref.repeat.bed"
    output: 
        path_names["pamref"]+"/unique_repeats.tsv"
    shell:
        "cat {input} | cut -f 4 | sed 's/::/\\t/g' | cut -f 2 |  sort | uniq -c > {output}"


rule format_repeat_masking:
    input:
        "{sample}.fa.out"
    output:
        "{sample}.repeat.bed"
    shell:
        "cat {input} | awk 'BEGIN{{OFS=\"\\t\";}}NR>3{{gsub(\"C\",\"-\",$9);print $5,$6-1,$7,$10\"::\"$11,0,$9;}}' > {output}"

rule make_repeat_mask_gff:
    input:
        path_names["pamref"]+"/events_ref.fa"
        #config["results"]+"/events.fa" 
    output:
        path_names["pamref"]+"/events_ref.fa.out"
        #config["results"]+"/events.fa.out" 
    threads:
        config["repeatmasker_threads"]
    shell:
        "RepeatMasker -species human -pa {threads} {input}"



"""
{
    "# of Primary Mappings": "1783597",
    "# of Supp. Mappings": "2",
    "Chimeric": "0",
    "Concordant": "1749040",
    "Discordant": "609",
    "Mean": "499.50",
    "OEA": "9509",
    "Orphan": "24439",
    "Range": "[454, 545]",
    "Read Length Range": "[100, 100]",
    "Std": "26.56",
    "Total Number of Reads": "1783597"
}

"""
rule make_stats_js:
    input:
        expand(path_names["remcor"] + "/{sample}/{sample}.json",sample=config["input"].keys()),
    output:
        config["results"]+"/stats.js"
    run:
        import json
        import sys
        buff= []
        short = {
            "# of Primary Mappings": "p",
            "# of Supp. Mappings": "s",
            "Chimeric": "c",
            "Concordant": "C",
            "Discordant": "D",
            "Mean": "m",
            "OEA": "E",
            "Orphan": "O",
            "Range": "r",
            "Read Length Range": "l",
            "Std": "v",
            "Total Number of Reads": "t"
        }

        for f in input:
            with open( f , 'r') as hand:
                end = f.rfind(".json")
                start = f.rfind("/") + 1
                sample = f[start:end]
                jsn = json.load(hand)
                parts = []
                for k,v in jsn.items():
                    if k not in short:
                        continue
                    parts.append('"{}":"{}"'.format(short[k],v))
                
                buff.append('"{}":{}'.format(sample,"{"+ ",".join(parts) + "}"))
        with open(output[0], 'w') as hand:
            print("/*{}*/".format(",".join([":".join([k,v]) for k,v in short.items()])),file=hand) # Short to proper comment
            print("stats = ",end="",file=hand) # variable
            print("{" + ",".join(buff) + "}",file=hand) # fields


rule make_data_js:
    input:
        expand(path_names["pamgen"]+"/{sample}/final.vcf",sample=config["input"].keys()),
    output:
        config["results"]+"/data.js"
    run:
        import json
        import sys
        populations = {}
        cohort = config["population"]
        sn = {"uuid" :"id","qual" :"q","genotype" :"g","Cluster" :"c","Support" :"p","FLSUP" :"flp","FRSUP" :"frp","FSUP" :"fsp","GLeft" :"L","GRight" :"T","GRef" :"R","GLRatio" :"glr","GRRatio" :"gtr","Sample" :"s", "unique_content": "u", "repeat_ranges": "r" }
        unwanted_tag_list=["FLSUP","FRSUP","FSUP","GLRatio","GRRatio", "repeat_ranges", "unique_content"] 
        with open(output[0], 'w') as ohand:
            print("/*\nid uuid\nq qual\ng genotype\nc Cluster\np Support\nL GLeft\nT GRight\nR GRef\ns Sample\n*/",file=ohand)
            print("var mut_data = ",file=ohand)
            for invcf in input:
                sampleid = invcf.split("/")[-2]
                if cohort not in populations:
                    populations[cohort] = {"individuals" : {}}
                
                if sampleid not in populations[cohort]["individuals"]:
                    populations[cohort]["individuals"][sampleid] = {'events' : None}
                variants = []
                with open(invcf,'r') as ihand:
                    for line in ihand:
                        if line[0] == "#":
                            continue
                        fields = line.rstrip().split("\t")
                        variant = {}

                        variant["id"] = fields[2]        
                        variant["q"] = "{:.1f}".format(float(fields[5]))  
                        variant["g"] = fields[10]        
                        info = fields[7].split(";")
                        for i in info:
                            parts = i.split("=")
                            if len(parts) == 2:
                                if parts[0] in unwanted_tag_list:
                                    continue
                                else:
                                    variant[sn[parts[0]]] = parts[1]
                            else:
                                variant[parts[0]] = ""

                        variants.append(variant)
                populations[cohort]["individuals"][sampleid]["events"] = variants
                
            ohand.write(json.dumps(populations))
            ohand.write(";")

rule all_rem_temp_vcfs:
    input:
        vcf=expand(path_names["pamass"] +"/{sample}/all.vcf",sample=config["input"].keys())
    output:
        temp(config["results"] + "/rem.done"),
    shell:
        "rm -f " + path_names["pamass"] + "/*/T*; touch {output}"
   

rule make_summary_js:
    input:
        table=path_names["pamgen"]+"/table.tsv",
    output:
        config["results"]+"/summary.js"
    run:
        def sumall(l):
            sm = 0
            for i in l:
                sm+=i
            return sm

        import sys
        import json
        length2counts = {}
        repeats = {}
        with open(output[0],'w') as ohand, open(input.table,'r') as ihand:
            print("/*event\te\ncount\tc\nlength\tl\nrepeat_type\tr\nunique_content_percentage\tu\n*/",file=ohand)
            print("var cohort_request=\"{}\"".format(config["population"]),file=ohand)
            print("var cohort_size={}".format(len(config["input"].keys())),file=ohand)
            print("var table_summary_file =[",file=ohand)
            for line in ihand:
                fields=line.rstrip().split("\t")
                print(json.dumps({"e" : fields[1], "c" : str(int(float(fields[0]))), "l":len(fields[2]), "r":fields[3], "u":fields[4]}),file=ohand,end=",")
                if fields[3] not in repeats:
                    repeats[fields[3]] = 1
                length = int(len(fields[2])/5)*5
                count = int(float(fields[0]))
                if length not in length2counts:
                    length2counts[length] = []
                length2counts[length].append(count)
            print("];",file=ohand)
            sorted_lens = sorted(length2counts.keys())
            print("var length_labels = [",file=ohand)
            print( ", ".join([str(x) for x in sorted_lens]),file=ohand,end=",")
            print("];",file=ohand)


            print("var unique_counts = [",file=ohand)
            print( ", ".join([ str(len(length2counts[x])) for x in sorted_lens]),file=ohand,end=",")
            print("];",file=ohand)

            print("var total_counts = [",file=ohand)
            print( ", ".join([ str(sumall(length2counts[x])) for x in sorted_lens]),file=ohand,end=",")
            print("];",file=ohand)
            

            print("var repeat_types = [",file=ohand)
            reps = []
            for line in repeats.keys():
                reps.append("\"{}\"".format(line))
            print(", ".join(reps),end="",file=ohand)
            print("];",file=ohand)

rule annotate_repeats_in_vcf:
    input:
        bed=path_names["pamref"]+"/ind/{sample}/events_ins_pos_and_repeat.bed",
        vcf=path_names["pamgen"]+"/{sample}/genotyped.vcf",
    output:
        vcf=path_names["pamgen"]+"/{sample}/final.vcf",
    run:

        bedl = []
        vcfl = []
        with open(input.vcf , 'r') as vhand, open(input.bed, 'r') as bhand, open(output.vcf, 'w') as ohand:
            bedln = bhand.readline()
            vcfln = vhand.readline()
            while(bedln and vcfln):
                while(vcfln[0] == "#"):
                    print(vcfln,end="",file=ohand)
                    vcfln = vhand.readline()
               
                vcffields = vcfln.rstrip().split("\t")
                bedfields = bedln.rstrip().split("\t")
                bedl.append(bedfields)
                vcfl.append(vcffields)
                bedln = bhand.readline()
                vcfln = vhand.readline()
                
            bedl = sorted(bedl, key=lambda x: x[0])
            vcfl = sorted(vcfl, key=lambda x: x[2])

            for b,v in zip(bedl,vcfl):
                print("\t".join(v[:8]),end=";",file=ohand)
                bb = b[3].split(";")
                print(bb[-2],bb[-1],sep=";",file=ohand,end="\t") 
                print("\t".join(v[8:]),file=ohand)

rule all_vis_table:
    input:
        expand(path_names["pamgen"]+"/{sample}/final.vcf",sample=config["input"].keys()) 
    output:
        path_names["pamgen"]+"/table.tsv"
    run:
        event_details = {}
        events = {}
        total_count = len(input)
        for vcf in input:
            my_id = vcf.split("/")[-2]
            with open( vcf , 'r') as hand:
                for line in hand:
                    if line[0] == "#":
                        continue
                    fields = line.rstrip().split("\t")
                    infocol = fields[7].split(";")
                    info = {}
                    for fi in infocol:
                        parts = fi.split("=")
                        if len(parts) == 1:
                            info[parts[0]] = ""
                        else:
                            info[parts[0]] = parts[1]
                    #print(info)
                    event_id = fields[2]
                    if event_id not in events:
                        events[event_id] = []
                        event_details[event_id] = (fields[2],fields[4][1:],",".join([x.split("(")[0] for x in info["repeat_ranges"].split(":")]), info["unique_content"])
                    if fields[-1] != "0/0" and fields[6] == "PASS":
                        events[event_id].append((my_id,fields[-1]))
        with open( output[0], 'w') as hand:
            for event,patients in events.items():
                details = event_details[event]
                cnt = float(len(patients))
                if cnt < 1:
                    continue
                print( "{}\t{}".format(cnt,"\t".join(details)), file=hand)


rule get_final_bam:
    input:
        bam=path_names["pamgen"]+"/{sample}/{sample}.bam",
        orphan_bam=path_names["pamgen"]+"/{sample}/orphans.bam",
    output:
        orphan_bam=path_names["pamgen"]+"/{sample}/final.bam",
    threads:
        config["other_threads"]
    shell:
        "samtools merge {output} {input.bam} {input.orphan_bam} -c -@ {threads} && samtools index -@ {threads} {output}"

rule convert_vis_sam_to_bam:
    input:
        sam=path_names["pamgen"]+"/{sample}/{sample}.sam",
        header=path_names["pamgen"]+"/{sample}/header.sam",
    output:
        temp(path_names["pamgen"]+"/{sample}/{sample}.bam"),
    threads:
        16
    shell:
        "cat {input.header} {input.sam} | samtools sort -m6G -@ {threads} > {output}"

rule make_rg_lines_for_header:
    output:
        temp(path_names["pamgen"]+"/rg.head"),
    run:
        with open(output[0], 'w') as hand:
            categories = ["BreakPointSupport",
                    "BreakPointSupportOEA",
                    "BreakPointSupportOrphan",
                    "ReferenceSupport",
                    "FalseSupport",
                    "NoSupport"]
            categories = categories + [ "DIS-{}".format(x) for x in categories]
            for cat in categories:
                print("@RG\tID:{0}\tLB:{0}".format(cat),file=hand)

rule genotype_vis:
    input:
        cram=get_cram_name,
        cram_index=get_cram_index,
        vcf=path_names["pamscv"]+"/{sample}/setcover-filtered.named.no-centro.vcf",
        fasta=path_names["pamref"]+"/events_ref.fa",
        read_stats=path_names["remcor"] +"/{sample}/{sample}.stats.json",
    output:
        sam=temp(path_names["pamgen"]+"/{sample}/{sample}.sam"),
        vcf=temp(path_names["pamgen"]+"/{sample}/genotyped.vcf"),
        head=temp(path_names["pamgen"]+"/{sample}/{sample}.header"),
        bed  =temp(path_names["pamref"]+"/ind/{sample}/events_ins_pos.bed"),
    params:
        ref=config["reference"],
        rang=config["genotyping_flank_size"],
    threads:
        config["other_threads"]
    run:

        import json
        with open(input.read_stats, 'r') as jhand:
            read_stats =json.load(jhand)

        rl = read_stats["read_len"]
        min_frag = read_stats["tlen_min"]
        max_frag = read_stats["tlen_max"]
 
        genotype_dic = { 0 : "0/0", 1 : "1/0", 2 : "1/1"}

        cmd="samtools view -H {} > {}".format(input.cram,output.head)
        process = subprocess.Popen(cmd,shell=True);
        process.communicate()
        view_cmd="samtools view -T {} {}".format(params.ref,input.cram)
        sort_cmd = "samtools sort -n | samtools view"
        in_house_cmd = SCPA + "/pamir process_range {} {} {} {} {} {} {}"

        id2seq = {}

        with open(input.fasta,"r") as hand:
            line = hand.readline()
            while line:
                fields = line.rstrip().split()
                name = fields[0][1:]
                line = hand.readline()
                id2seq[name] = (line.rstrip(),fields[1],fields[2])
                line = hand.readline()


        with open(input.vcf,"r") as hand , open( output.sam, "w") as samhand, open(output.vcf, "w") as vcfhand, open(output.bed, "w") as bedhand:



            line = hand.readline()
            while line:
                while line[0] == "#":
                    print(line.rstrip(),file=vcfhand)
                    line = hand.readline()

                fields = line.rstrip().split("\t")

                processes= []
                ps = []
                pm = [] 
                vcfs = []
                beds = []
                for tid in range(threads):
                    while line and fields[6] != "PASS":
                        line = hand.readline()
                        fields = line.rstrip().split("\t")
                    if not line:
                        break
                    pos = int(fields[1])
                    start = pos - params.rang + 1
                    end   = pos + params.rang - 1
                    seq = fields[4][1:]
                    vcfs.append(fields)
            #fasta_cmd="cat  {} | grep {} -A1".format(input.fasta, fields[2])
            #ps.append( subprocess.Popen(fasta_cmd, shell=True, stdout=subprocess.PIPE))
                    ps.append(fields[2])
                    pm.append(view_cmd + " {}:{}-{}".format(fields[0],start+rl,end-rl) + "|" + in_house_cmd.format(pos,len(seq),params.rang,min_frag,max_frag,fields[2],"{}"))
                    line = hand.readline()
                    fields = line.rstrip().split("\t")
                pf = []
                for name,m in zip(ps,pm):

                    if name not in id2seq:
                        continue
                    tup = id2seq[name]

                    fasta = tup[0]
                    bedinfo = tup[1:]


                    beds.append("{}\t{}\t{}".format(name,bedinfo[0],int(bedinfo[0]) + int(bedinfo[1])))

                    pf.append(m.format(fasta))

                    

                for f in pf:
                    processes.append(subprocess.Popen(f, shell=True,stdout=subprocess.PIPE))


                for bed,v,process in zip(beds,vcfs,processes):
                    stdout,stderr = process.communicate()
                    samlines = stdout.decode('utf-8').splitlines()

                    for l in samlines[1:]:
                        print(l,file=samhand)
                    print("\t".join(v[:8]),end=";",file=vcfhand)

                    #print(samlines[0])
                    left,right,ref = [int(x) for x in samlines[0].split("\t")]
                    
                    if ref + left == 0:
                        left_ratio  = -1
                    else:
                        left_ratio  = (left  - ref) / (left  + ref)
                    if ref  + right == 0:
                        right_ratio = -1
                    else:
                        right_ratio = (right - ref) / (right + ref)

                    lg = []
                    for ratio in [left_ratio, right_ratio]:
                        if ratio >= 0.3:
                            lg.append(2)
                        elif ratio <= -0.3:
                            lg.append(0)
                        else:
                            lg.append(1)
                    
                    genotype = genotype_dic[min(lg)]

                    infoformat="GLeft={};GRight={};GRef={};GLRatio={:.2f};GRRatio={:.2f};Sample={}"

                    print( infoformat.format(left,right,ref,left_ratio,right_ratio,wildcards.sample),end="\t",file=vcfhand)
                    print( bed,end="\t",file=bedhand)
                    print( infoformat.format(left,right,ref,left_ratio,right_ratio,wildcards.sample),end="\n",file=bedhand)

                    print("\t".join(v[8:]),file=vcfhand,end="\t")
                    print("GT\t{}".format(genotype),file=vcfhand)

                if not line:
                    break

rule get_insertions_sam_header:
    input:
        bam=path_names["pamgen"]+"/{sample}/orphans.bam",
    output:
        temp(path_names["pamgen"]+"/{sample}/header.sam"),
    shell:
        "samtools view {input} -H > {output}"
       
rule convert_orphans_to_bam:
    input:
        sam=path_names["pamgen"]+"/{sample}/orphans.sam",
        rgs=path_names["pamgen"]+"/rg.head",
    output:
        bam=temp(path_names["pamgen"]+"/{sample}/orphans.bam"),
    threads:
        16
    shell:
        "cat <(samtools view -H {input.sam}) {input.rgs} <(samtools view {input.sam}) | samtools view -Shb | samtools sort -m4G -@ {threads} -o {output.bam}"
rule map_orphans_to_vis_fasta:
    input:
        fasta=path_names["pamref"]+"/events_ref.fa",
        fasta_index=path_names["pamref"]+"/events_ref.fa.bwt",
        fastq=path_names["remcor"] + "/{sample}/{sample}.orphan.canonical."+sext,
    output: 
        sam=temp(path_names["pamgen"]+"/{sample}/orphans.sam"),
    threads:
        config["aligner_threads"]
    params:
        scpa=SCPA
    shell:
        "bwa mem -t {threads} {input.fasta} {input.fastq} -p | {params.scpa}/pamir process_orphan > {output}"

rule bwa_index:
    input:
        fasta="{any}.fa",
    output:
        index="{any}.fa.bwt",
    shell:
        "bwa index {input}"


rule merge_fastas:
    input:
        expand(path_names["pamref"]+"/ind/{sample}/events_ref.fa",sample=config["input"].keys())
    output:
        fasta=path_names["pamref"]+"/events_ref.fa",
        fasta_index=path_names["pamref"]+"/events_ref.fa.fai",
        fasta_lookup=path_names["pamref"]+"/events_ref.fa.lookup",
    params:
        scpa=SCPA,
    shell:
        "printf \"{input}\n\" | tr ' ' '\n' | python {params.scpa}/scripts/merge_refs.py {output.fasta_lookup} > {output.fasta} && samtools faidx {output.fasta}"
   
rule make_vis_fasta:
    input:
        vcf=path_names["pamscv"]+"/{sample}/setcover-filtered.named.no-centro.vcf",
    output:
        fasta=temp(path_names["pamref"]+"/ind/{sample}/events_ref.fa"),
    params:
        ref=config["reference"],
        flank=config["genotyping_flank_size"],
    run:
        import subprocess
        contig_sizes = {}
        with open("{}.fai".format(params.ref), 'r') as fasta_index:
            for line in fasta_index:
                fields = line.split("\t")
                contig_sizes[fields[0]] = int(fields[1])
        with open(input['vcf'],'r') as hand, open(output["fasta"],"w+") as fast_hand:
            for line in hand:
                if line[0] == "#":
                    continue
                fields = line.rstrip().split("\t")
                if fields[6] != "PASS":
                    continue
 #               if fields[9] == "0/0":
 #                   continue
                ins_id = fields[2]
                ch = fields[0]
                pos = int(fields[1])
                mid = params.flank
                start = pos - params.flank
                end   = pos + params.flank 
           
                if start < 0:
                    mid = mid + start
                    start = 0
                if end > contig_sizes[ch]:
                    end = contig_sizes[ch] - 1
                cmd = "printf \"{}\\t{}\\t{}\\n\" | bedtools getfasta -fi {} -bed stdin".format(ch,start,end,params["ref"])

                ins = fields[4][1:]
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                stdout,stderr = process.communicate()
                stdout = stdout.decode('utf-8').splitlines()
                
                #print("{}\t{}".format(cmd,fields[2]))
                seq = stdout[-1]
                seq = seq[:mid] + ins + seq[mid:]
                print(">{} {} {}\n{}".format(ins_id, mid, len(ins),seq.rstrip()),file=fast_hand)


rule all_filtered_vcf:
    input:
        expand(path_names["pamscv"]+"/{sample}/setcover-filtered.vcf",sample=config["input"].keys()) 
    output:
        path_names["pamann"]+"/done"
    shell:
        "touch {output}"

if "centromeres" in config:
    rule remove_centromeres_from_vcf:
        input:
            "{sample}.vcf",
        output:
            "{sample}.no-centro.vcf",
        params:
            centromeres=config["centromeres"],
        shell:
            "bedtools intersect -a {input} -b <(cat {params.centromeres} | sed 's/chr//') -v -header > {output}"
else:
    rule mock_remove_centromeres_from_vcf:
        input:
            "{sample}.vcf",
        output:
            "{sample}.no-centro.vcf",
        shell:
            "cp {input} {output}"

rule rename_events:
    input:
        "{sample}.vcf"
    output: 
        temp("{sample}.named.vcf")
    run:
        import hashlib
        with open(input[0], 'r') as ihand, open(output[0], 'w') as ohand:
            for line in ihand:
                if( line[0] == "#"):
                    print(line.rstrip(), file=ohand)
                    continue
                fields = line.rstrip().split("\t")
                ho = hashlib.md5(fields[4].encode())
                fields[2] = "-".join([fields[0],fields[1],ho.hexdigest()[:11]])
                print( "\t".join(fields), file=ohand)
               
rule filter_by_setcover:
    input:
        smooth=path_names["pamscv"]+"/{sample}/smooth",
        fvcf=path_names["pamann"]+"/{sample}/filtered.vcf",
    output:
        temp(path_names["pamscv"]+"/{sample}/setcover-filtered.vcf"),
    params:
        scpa=SCPA,
    shell:
        "python {params.scpa}/scripts/filter_by_setcover.py {input.smooth} {input.fvcf}  {output}"

rule smoother:
    input:
        path_names["pamscv"]+"/{sample}/prep",
    output:
        temp(path_names["pamscv"]+"/{sample}/smooth"),
    params:
        scpa=SCPA
    shell:
        "{params.scpa}/pamir smoother {input} >  {output}"


if new_prep:
    rule prep_for_set_cover:
        input:
            fsc=path_names["pamann"] +"/{sample}/filtered.vcf",
            sam=path_names["pamann"]+"/{sample}/seq.mrsfast.recalsorted.sam",
        output:
            temp(path_names["pamscv"]+"/{sample}/prep"),
        params:
            scpa=SCPA,
        shell:
            "python {params.scpa}/scripts/new_generate_setcover.py {input.sam} {input.fsc} {output}"
else: 
    rule prep_for_set_cover:
        input:
            fsc=path_names["pamann"] +"/{sample}/all.sorted.vcf_filtered_for_setcov.sorted",
            sam=path_names["pamann"]+"/{sample}/seq.mrsfast.recalsorted.sam",
        output:
            temp(path_names["pamscv"]+"/{sample}/prep"),
        params:
            scpa=SCPA,
        shell:
            "python {params.scpa}/scripts/generate_setcover_input.py {input.fsc} {input.sam} {output}"
 
rule filter_vcf:
    input:
        vcf=path_names["pamass"] +"/{sample}/all.sorted.vcf",
        fq=path_names["remcor"]+"/{sample}/{sample}.all_interleaved."+sext,
        header=path_names["pamann"]+"/header.vcf",
        statfile=path_names["remcor"]+"/{sample}/{sample}.stats.json",
    output:
        vcf=path_names["pamann"]+"/{sample}/filtered.vcf",
        vcf_nh=temp(path_names["pamann"]+"/{sample}/filtered.vcf_nohead"),
        sam=path_names["pamann"]+"/{sample}/seq.mrsfast.recalsorted.sam",
    params:
        ref=config["reference"],
        tlen=1000,
        wd=path_names["pamann"],
        scpa=SCPA,
    threads:
        config["aligner_threads"]
    shell:
        "python {params.scpa}/scripts/filtering.py {input.vcf} {params.ref} {params.wd}/{wildcards.sample}/ {params.tlen} {threads} {input.fq} {input.statfile}  {output.vcf_nh} {output.sam} && cat {input.header} {output.vcf_nh} > {output.vcf}"
         
rule sort_vcf:
    input:    
        "{sample}.vcf"
    output:
        "{sample}.sorted.vcf"
    params:
        scpa=SCPA,
        dupc=1,
    shell:
        "python {params.scpa}/scripts/sort_vcf.py {input} {output} {params.dupc}"



    
rule generate_vcf_header:
    output:
        path_names["pamann"]+"/header.vcf",
    params:
        ref=config["reference"],
    run:
        header_info = "".join(["##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##reference={}\n##source=Pamir\n".format(params.ref),
                        "##INFO=<ID=Cluster,Number=1,Type=Integer,Description=\"ID of the cluster the variant is extracted from\">\n",
                        "##INFO=<ID=Support,Number=1,Type=Integer,Description=\"Number of reads/contigs supporting the contig\">\n",
                        "##INFO=<ID=SEQ,Number=1,Type=String,Description=\"Variant sequence\">\n",
                        "##INFO=<ID=FLSUP,Number=1,Type=Integer,Description=\"Number of left supporting reads in filtering\">\n",
                        "##INFO=<ID=FLRSUP,Number=1,Type=Integer,Description=\"Number of right supporting reads in filtering\">\n",
                        "##INFO=<ID=FSUP,Number=1,Type=Integer,Description=\"Number of total supporting reads in filtering\">\n",       
                        "##INFO=<ID=GLeft,Number=1,Type=Integer,Description=\"Number of reads passing through the left breakpoint on template INS sequence\">\n",
                        "##INFO=<ID=GRight,Number=1,Type=Integer,Description=\"Number of reads passing through the right breakpoint on template INS sequence)\">\n",
                        "##INFO=<ID=GRef,Number=1,Type=Integer,Description=\"Number of reads passing through the breakpoint on template REF sequence)\">\n",
                        "##INFO=<ID=GRRatio,Number=1,Type=Float,Description=\"Ratio between GRight and GRef\">\n",
                        "##INFO=<ID=GLRatio,Number=1,Type=Float,Description=\"Ratio between GLeft and GRef \">\n",
                        "##INFO=<ID=RMContent,Number=1,Type=String,Description=\"Repeat content overlapping with the eveRepeat content overlapping with the event\">\n",##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGENOTYPE\n"])
        with open(output[0],'w') as hand:
            print(header_info, file=hand, end="")
    

   
rule all_vcf:
    input:
        expand(path_names["pamass"] +"/{sample}/all.vcf",sample=config["input"].keys()) 
    output:
        path_names["pamass"] +"/merged.vcf"
    shell:
        "cat {input} > {output}"

rule all_cc:
    input:
        expand(path_names["partition"] +"/{sample}/cc",sample=config["input"].keys())
    output:
        path_names["partition"] +"/done"
    shell:
        "touch {output}"


#print(sext)
def get_assembly_input():
    return expand( "{Path}/{sample}/{sample}.orphan.{typ}."+sext, Path=path_names["remcor"], sample=config["input"].keys(),typ=["canonical","almost"])

def format_spades_inputs():
    return " ".join([ " --pe{0}-12 {2}/{3}/{3}.orphan.canonical.{4} --pe{1}-12 {2}/{3}/{3}.orphan.almost.{4}".format(1+i*2,1+i*2+2,path_names["remcor"],x,sext) for i,x in enumerate(config["input"].keys())])

if config["assembler"] == "minia":
    rule minia_cf:
        input:    
            get_assembly_input()
        output:
            path_names["extass"] + "/reads.fofn",
        params:
            analysis=config["analysis"],
        shell:
            "ls {}/*/*.orphan.*.{} > {{output}}".format(path_names["remcor"],sext)

    rule minia_all:
        input:     
            path_names["extass"] + "/reads.fofn",
        output:
            path_names["extass"] + "/contigs.fasta"
        params:
            k=config["assembler_k"],
            min_abundance=config["minia_min_abundance"],
            max_memory=250000,
            dr=path_names["extass"] + "/",
        threads:
            config["assembly_threads"]
        shell:
           "cd {params.dr} && minia  -verbose 0 -in {input} -kmer-size {params.k} -abundance-min {params.min_abundance} -nb-cores {threads} && mv reads.contigs.fa contigs.fasta && rm -f reads.unitigs.fa* reads.h5"
elif config["assembler"] == "spades": 
    rule spades_all:
        input:     
            get_assembly_input()
        output:
            path_names["extass"] + "/contigs.fasta"
        params:
            k=config["assembler_k"],
            min_abundance=10,
            max_memory=500,
            dr=path_names["extass"],
            formatted_input=format_spades_inputs(),
        threads:
            config["assembly_threads"]
        shell:
            "spades.py -m {params.max_memory} -t {threads} {params.formatted_input} -o {params.dr} && cat {params.dr}/contigs.fasta | tr  '\\n' '\\t' | sed 's/>/\\n>/g' | sed 's/\\t/\\n/' | tr -d '\\t' | awk 'NR>1' > {params.dr}/tmp.fa && mv {params.dr}/tmp.fa {output}"
elif config["assembler"] == "abyss":
     rule abyss_all:
        input:
            expand(path_names["remcor"]+"/{sample}/{sample}.orphan.{orphan_type}." + sext ,orphan_type=["canonical", "almost"],sample=config["input"].keys()),
            #expand(path_names["remcor"]+"/{sample}/{sample}.orphan.canonical." + sext + " " + path_names["remcor"] + "/{sample}/{sample}.orphan.almost." + sext,sample=[x[0][0:x[0].find(".")] for x in config["input"].values()]),
        output:
            path_names["extass"] + "/contigs.fasta"
        params:
            k=config["assembler_k"],
            dr=path_names["extass"] + "/",
        threads:
            config["assembly_threads"]
        shell:
            "cd {params.dr} && abyss-pe name=temp in='{input}' k={params.k} && mv temp-contigs.fa contigs.fasta"
elif config["assembler"] == "present":
    assert "assembly" in config, '"assembly" should be present in config to use present option'
    rule present_assembly_copy:
        input:
            config["assembly"]
        output:
            path_names["extass"] + "/contigs.fasta"
        shell:
            "cp {input} {output}"

rule pamir_assemble_full_new:
    input:
        partition=path_names["partition"] +"/{sample}/partition",
        cluster_count=path_names["partition"] +"/{sample}/partition.count",
        read_stats=path_names["remcor"] + "/{sample}/{sample}.stats.json",
    output:
        vcf=temp(path_names["pamass"] +"/{sample}/all.vcf"),
        vcf_lq=temp(path_names["pamass"] +"/{sample}/all_LOW_QUAL.vcf"),
        logs=path_names["pamass"] +"/{sample}/all.log",
    params:
        pppt = config["pamir_partition_per_thread"],
        ref=config["reference"],
        wd=path_names["pamass"],
    threads:
        config["assembly_threads"],
    run:
        with open( input.cluster_count, 'r') as chand:
            cc = int(chand.readline())
        
        index = 1
        use_threads = max(1,threads - 1)
        pppt = min(params.pppt,int(cc/use_threads))
       
        import json
        with open(input.read_stats, 'r') as jhand:
            read_stats =json.load(jhand)
        rl = read_stats["read_len"]
        min_frag = read_stats["tlen_min"]
        max_frag = read_stats["tlen_max"]

        

        cmd_template =  SCPA + "/pamir assemble {0} {1} {{0}}-{{1}} T{{2}} 30000 {2} {3}/{4} > {3}/{4}/T{{2}}.log".format(input.partition, params.ref, rl, params.wd, wildcards.sample)
        while index + pppt -1 <= cc:
            tids = []
            procs = []
            first_index = index
            for tid in range(use_threads):
                start = index
                if index + pppt > cc:
                    end = cc
                else:
                    end = index + pppt - 1
                cmd = cmd_template.format(start,end,tid)
                procs.append(subprocess.Popen(cmd, shell=True))
                tids.append(tid)
                index = end + 1

            print("Processing partitions between {} and {} with {} threads".format(first_index,index,len(procs),file=sys.stderr))
            for i,proc in enumerate(procs):
                proc.communicate()
            vcfs = [ "{}/{}/T{}.vcf".format(params.wd,wildcards.sample,tid) for tid in tids]
            lqvs = [ "{}/{}/T{}_LOW_QUAL.vcf".format(params.wd,wildcards.sample,tid) for tid in tids]
            logs = [ "{}/{}/T{}.log".format(params.wd,wildcards.sample,tid) for tid in tids]
            

            catcmd = "cat " + " ".join(vcfs) + " >> " + output.vcf + ".tmp"
            vcf_p  = subprocess.Popen(catcmd,shell=True)
            catcmd = "cat " + " ".join(lqvs) + " >> " + output.vcf_lq + ".tmp"
            low_p  = subprocess.Popen(catcmd,shell=True)
            catcmd = "cat " + " ".join(logs) + " >> " + output.logs + ".tmp"
            log_p  = subprocess.Popen(catcmd,shell=True)

            vcf_p.communicate()
            low_p.communicate()
            log_p.communicate()
  
        mvcmd = "mv " +  output.vcf + ".tmp " + output.vcf
        vcf_p  = subprocess.Popen(mvcmd,shell=True)

        mvcmd = "mv " +  output.vcf_lq + ".tmp " + output.vcf_lq
        low_p  = subprocess.Popen(mvcmd,shell=True)

        mvcmd = "mv " +  output.logs + ".tmp " + output.logs
        log_p  = subprocess.Popen(mvcmd,shell=True)

        vcf_p.communicate()
        low_p.communicate()
        log_p.communicate()
  

rule recalibrate_oea_to_orphan:
    input:
        lookup=path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.merged.fa.lookup",
        bam=path_names["pamoea"] + "/{sample}/{sample}.oea2orphan.bam"
    output:
        path_names["pamoea"] + "/{sample}/{sample}.oea2orphan.recalib.sam"
    shell:
        SCPA + "/pamir recalibrate {input.lookup} <(samtools view {input.bam}) {output}"

rule pamir_partition:
    input:
        sam=path_names["pamoea"] + "/{sample}/{sample}.anchor.sorted.sam",
        oea2orphan=path_names["pamoea"] + "/{sample}/{sample}.oea2orphan.recalib.sam",
        oea_unmapped=path_names["remcor"] + "/{sample}/{sample}.oea.unmapped."+sext,
        contigs=path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.fa",
    output:
        partition=path_names["partition"] +"/{sample}/partition",
        cluster_count=path_names["partition"] +"/{sample}/partition.count",
    params:
        rang=1000,
    shell:
        SCPA + "/pamir partition {input.sam} {output.partition} {params.rang} {input.contigs} {input.oea2orphan} {input.oea_unmapped}"


rule merge_contigs:
    input:
        path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.fa"
    output:
#         single=config["analysis"]+"/"+config["assembler"]+"/{sample}.reads.contigs.filtered.clean.single.fa",
        merged=path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.merged.fa",
        index=path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.merged.fa.lookup",
    params:
        scpa=SCPA,
    shell:
        "python {params.scpa}/scripts/prep-ctgs.py {input} {output.merged} 10"

if "blastdb" in config:
    rule blast_contigs:
        input:
            path_names["extass"] + "/contigs.fasta"
        output:
            path_names["blast"] + "/contigs.megablast"
        threads:
            config["other_threads"],
        params:
            db=config["blastdb"],
        shell:    
            "blastn -task megablast -db {params.db} -num_threads {threads} -query {input} -outfmt 6 > {output}"

    rule find_contaminations:
        input:
            path_names["blast"] + "/contigs.megablast",
        output:
            path_names["blast"] + "/contigs.megablast.tabular",
        shell:
            "cat {input} | cut -f 1 > {output}"

    rule clean_contigs:
        input:
            fasta=path_names["extass"] + "/contigs.fasta",
            mb=path_names["blast"] + "/contigs.megablast.tabular",
        output:
            path_names["blast"] + "/contigs.clean.fa"
        params:
            scpa=SCPA,
        shell:
            "python {params.scpa}/scripts/remove_contaminations.py {input.mb} {input.fasta} {output}"

else:
    rule mock_clean_contigs:
        input:
            fasta=path_names["extass"] + "/contigs.fasta",
        output:
            path_names["blast"] + "/contigs.clean.fa"
        shell:
            "ln -s {input} {output}"

rule filter_contigs:
    input:
        fa=path_names["blast"] + "/contigs.clean.fa",
        json=path_names["remcor"] + "/{sample}/{sample}.stats.json",
    output:
        temp(path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.fa")
    run:
        if "pamir_min_contig_len" in config:
            min_ctg_len = config["pamir_min_contig_len"]
        else:
            import json
            with open(input.json, 'r') as jhand:
                stat_cfg = json.load(jhand)
            min_ctg_len = stat_cfg["read_len"]

        with open(input.fa,'r') as hand, open(output[0],'w+') as arm:
            line = hand.readline()
            while line != '':
                if line[0] != '>':
                    raise Exception('multi line fastas not good')
                header = line.rstrip()
                line = hand.readline()
                if len(line) > min_ctg_len + 1: #Plus 1 for endline
                    print(header,file=arm)
                    print(line.rstrip(),file=arm)
                line = hand.readline()

rule size_filter_for_vcf:
    input:
        "{sample}.vcf"
    output:
        "{sample}.l{mylen}.vcf"
    run:
        with open(input[0],'r') as hand, open(output[0],'w') as out:
            for line in hand:
                if line[0] == "#":
                    print(line.rstrip(),file=out)
                fields = line.rstrip().split("\t")
                info = fields[7].split(";")
                for i in info:
                    if "SVLEN" in i:
                        a,b = i.split("=")
                        if a > wildcards.mylen:
                            print(line.rstrip())
                            break

rule mrsfast_index_contigs:
    input:
        "{sample}.fa"
    output:
        "{sample}.fa.index"
    params:
        k=14,
    shell:
        "mrsfast --index {input} --ws {params.k}"

rule merge_oea_bams:
    input:
        full  = path_names["pamoea"] + "/{sample}/{sample}.oea.unmapped.sorted.bam",
        front = path_names["pamoea"] + "/{sample}/{sample}.oea.unmapped.front.sorted.bam",
        tail  = path_names["pamoea"] + "/{sample}/{sample}.oea.unmapped.tail.sorted.bam",
    output:
       temp(path_names["pamoea"] + "/{sample}/{sample}.oea2orphan.bam"),
    shell:
        "samtools merge {output} {input.full} {input.front} {input.tail} -@ {threads}" 

rule sam_sort:
    input:
        "{sample}.sam",
    output:
        "{sample}.sorted.sam"
    shell:
        "samtools sort  {input} -m 8G -@ {threads} -O SAM -o {output}"


rule bam_index:
    input:
        "{sample}.bam",
    output:
        "{sample}.bam.bai",
    shell:
        "samtools index {input}"


rule cram_index:
    input:
        "{sample}.cram",
    output:
        "{sample}.cram.crai",
    shell:
        "samtools index {input}"

rule bam_sort:
    input:
        "{sample}.sam",
    output:
        temp("{sample}.sorted.bam")
    shell:
        "samtools sort  {input} -m 8G -@ {threads} -o {output}"
 
rule mrsfast_oea_unmapped_contig_map_tail_cropped:
    input:
        nohit=path_names["pamoea"] + "/{sample}/{sample}.oea.unmapped.nohit.fq.gz",
	contigs=path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.merged.fa",
        index=path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.merged.fa.index",
    output:
        sam=temp(path_names["pamoea"] + "/{sample}/{sample}.oea.unmapped.tail.sam"),
    params:
        mem="8G",
        error=0,
        crop=75,
        N=50,
    threads:
        8
    shell:
        "mrsfast --disable-nohits --search {input.contigs} --seq {input.nohit} --threads {threads} -e {params.error} -o {output.sam} --mem {params.mem} --tail-crop {params.crop} -n {params.N}"+ mrsfast_comp_param 

def remove_gz(st):
	if st[-2:] == ".gz":
		return st[:-2]
	return st

rule mrsfast_oea_unmapped_contig_map_cropped:
    input:
        nohit=path_names["pamoea"] + "/{sample}/{sample}.oea.unmapped.nohit.fq.gz",
        contigs=path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.merged.fa",
        index=path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.merged.fa.index",
    output:
        sam=temp(path_names["pamoea"] + "/{sample}/{sample}.oea.unmapped.front.sam"),
    params:
        mem="8G",
        error=0,
        crop=75,
        N=50,
    threads:
        8
    shell:
        "mrsfast --disable-nohits --search {input.contigs} --seq {input.nohit} --threads {threads} -e {params.error} -o {output.sam} --mem {params.mem} --crop {params.crop} -n {params.N}"+ mrsfast_comp_param

rule mrsfast_oea_unmapped_contig_map:
    input:
        oea=path_names["remcor"] + "/{sample}/{sample}.oea.unmapped."+sext, 
        contigs=path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.merged.fa",
        index=path_names["lenfilter"] + "/{sample}.reads.contigs.filtered.clean.merged.fa.index",
    output:
        sam=temp(path_names["pamoea"] + "/{sample}/{sample}.oea.unmapped.sam"),
        nohit=path_names["pamoea"] + "/{sample}/{sample}.oea.unmapped.nohit.fq.gz"
    params:
        mem="8G",
        error=0,
	nh_path=lambda wildcards, output: remove_gz(output[1]),
    threads:
        8
    shell: 
        "mrsfast --search {input.contigs} --seq {input.oea}  --threads {threads} -e {params.error} -o {output.sam} -u {params.nh_path} --mem {params.mem}" + mrsfast_comp_param + " && gzip {params.nh_path}"

def check_reference_index(reference):
    import os.path
    index_exists = os.path.isfile(reference+".index")
    mrs_fast_index_cmd="mrsfast --index " + str(reference) + "; "
    return "echo Index file exists.; " if index_exists else mrs_fast_index_cmd

rule mrsfast_anchor_wg_map:
    input:
        path_names["remcor"] + "/{sample}/{sample}.oea.mapped."+sext, 
    output:
        temp(path_names["pamoea"] + "/{sample}/{sample}.anchor.sam"),
    params:
        mem="8G",
        ref=config["reference"],
        error=7,
        N=50,
        create_index=check_reference_index(config["reference"]),
    threads:
        config["aligner_threads"]
    shell:
        "{params.create_index} mrsfast --disable-nohits --search {params.ref} --seq {input} --threads {threads}  -e {params.error} -o {output} --mem {params.mem} -n {params.N}" + mrsfast_comp_param

rule velvet_all:
    input:
        expand(path_names["remcor"] + "/{sample}/{sample}.orphan.fq",sample=config["input"].keys())
    output:
        config["analysis"]+"/velvet/contigs.fasta"
    params:
        dr=config["analysis"]+"/velvet",
        velveth="velveth-omp",
        velvetg="velvetg-omp",
    threads:
        32
    shell:
        "module load velvet/1.2.10  && export NUM_OMP_THREADS {threads} && cd {params.dr}  && {params.velveth} . 31 -fastq <(cat {input}) && {params.velvetg} . -min_contig_lgth 400 -cov_cutoff 2; module unload velvet"

rule get_all_reads_from_cram:
    input:
        get_cram_name
    output:
        m1=config["analysis"]+"/fastqs/{sample}/{sample}.mate1.gz",
        m2=config["analysis"]+"/fastqs/{sample}/{sample}.mate2.gz",        
    params:
        ref=config["reference"],
        wd=path_names["remcor"] + "/{sample}",
        scpa = SCPA,
    shell:
        "cd {params.wd} && {params.scpa}/pamir getfastq <(samtools view -T {params.ref} {input}) {wildcards.sample} && mv {wildcards.sample}_1.fastq.gz {output.m1} && mv {wildcards.sample}_2.fastq.gz {output.m2}"

from os import path

def replace_sample_name(wildcards):
    return path.normpath(config["path"] + "/" + config["raw-data"] + "/" + config["input"][wildcards.sample][0])
rule link_bam:
    input:
        replace_sample_name
    output:
         path_names["linkeddata"]+"/{sample}.bam"
    shell:
        "ln -s {input} {output}"

rule link_cram:
    input:
        replace_sample_name
    output:
        path_names["linkeddata"]+"/{sample}.cram"
    shell:
        "ln -s {input} {output}"

"""
Total Number of Reads: 1790291
# of Primary Mappings: 1790291
  # of Supp. Mappings: 1

Original:
Concordant: 1748701
Discordant: 2396
  Chimeric: 143
       OEA: 10166
    Orphan: 28885

Processed:
Concordant: 1747025
Discordant: 806
  Chimeric: 0
       OEA: 11932
    Orphan: 30528

Read Length Range: [100, 100]

TLEN:
Range: [454, 545]
 Mean: 499.50
  Std: 26.56
"""

rule make_read_config:
    input:
        statfile=path_names["remcor"] + "/{sample}/{sample}.stat",
    output:
        read_stats=path_names["remcor"] + "/{sample}/{sample}.stats.json",
        seq_stats=path_names["remcor"] + "/{sample}/{sample}.json",
    run:
        import json
        stats_prep_cfg = {}       
        with open(input.statfile, 'r') as hand:
            for line in hand:
                fields = line.strip().split(": ")
                if len(fields) > 1:
                    stats_prep_cfg[fields[0]] = fields[1]


        rrange = stats_prep_cfg["Read Length Range"].replace("[","").replace("]","").split(", ")
        assert int(rrange[0]) == int(rrange[1]), "Error:\nVariable read length detected!\nRead lengths MUST be consistent within each input BAM/CRAM file!"
        stats_cfg = {}

        stats_cfg["read_len"] = int(rrange[0])
        if "tlen_min" in config:
            stats_cfg["tlen_min"] = config["tlen_min"]
        else:
            stats_cfg["tlen_min"] = int(float(stats_prep_cfg["Mean"]) - float(config["n_std_dev"]) * float(stats_prep_cfg["Std"]))  
        
        if "tlen_max" in config:
            stats_cfg["tlen_max"] = config["tlen_max"]
        else:
            stats_cfg["tlen_max"] = int(float(stats_prep_cfg["Mean"]) + float(config["n_std_dev"]) * float(stats_prep_cfg["Std"]))  

        with open(output.read_stats, 'w') as jhand:
            print(json.dumps(stats_cfg, indent=4, sort_keys=True), file=jhand)
        with open(output.seq_stats, 'w') as jhand:
            print(json.dumps(stats_prep_cfg,indent=4, sort_keys=True), file=jhand)

        
rule cram_split:
    input:
        get_cram_name
    output:
        orphan_can=path_names["remcor"] + "/{sample}/{sample}.orphan.canonical."+sext,
        orphan_ex=path_names["remcor"] + "/{sample}/{sample}.orphan.almost."+sext,
        oea_mapped=path_names["remcor"] + "/{sample}/{sample}.oea.mapped."+sext,
        oea_unmapp=path_names["remcor"] + "/{sample}/{sample}.oea.unmapped."+sext, 
        al=path_names["remcor"] + "/{sample}/{sample}.all_interleaved."+sext,
        statfile=path_names["remcor"] + "/{sample}/{sample}.stat",
    params:
        analysis=path_names["remcor"],
        REF=config["reference"],
        scpa=SCPA,
        PAMIR="./pamir",
        CRAMS= path_names["linkeddata"],
        pamir_params="{} 1 1 0.95".format(sext),
    threads:
        1
    shell:
        "CRPW={params.scpa} && cd {params.analysis}/{wildcards.sample} && samtools view {input} -T {params.REF} | $CRPW/{params.PAMIR} remove_concordant /dev/stdin {wildcards.sample} {params.pamir_params}"


