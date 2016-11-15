#!/usr/bin/env python
#786 

import os, sys, errno, argparse, subprocess, fnmatch, ConfigParser
#############################################################################################
# Class for colored texts and binary path
class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

class pipeline:
	pamir 	= os.path.dirname(os.path.realpath(__file__)) + "/pamir.py"
	sniper   	= os.path.dirname(os.path.realpath(__file__)) + "/sniper"
	sga		   	= os.path.dirname(os.path.realpath(__file__)) + "/sga.py"
	minia	   	= os.path.dirname(os.path.realpath(__file__)) + "/minia"
	velveth		= os.path.dirname(os.path.realpath(__file__)) + "/velveth"
	velvetg		= os.path.dirname(os.path.realpath(__file__)) + "/velvetg"
	mrsfast  	= os.path.dirname(os.path.realpath(__file__)) + "/mrsfast"
	bedtools  	= os.path.dirname(os.path.realpath(__file__)) + "/bedtools"
	samtools	= os.path.dirname(os.path.realpath(__file__)) + "/samtools"
	recalibrate = os.path.dirname(os.path.realpath(__file__)) + "/recalibrate"
	pprocessor  = os.path.dirname(os.path.realpath(__file__)) + "/partition_processor"
	sortvcf   = os.path.dirname(os.path.realpath(__file__)) + "/sort_vcf.py"
	ext_sup     = os.path.dirname(os.path.realpath(__file__)) + "/extract_support"
	filtering	= os.path.dirname(os.path.realpath(__file__)) + "/filtering.py"
	gensetcov	= os.path.dirname(os.path.realpath(__file__)) + "/generate_setcover_input.py"
	smoother	= os.path.dirname(os.path.realpath(__file__)) + "/smoother"
	genotyping	= os.path.dirname(os.path.realpath(__file__)) + "/allinone_genotyping.py"
	filterbysetcover	= os.path.dirname(os.path.realpath(__file__)) + "/filter_by_setcover.py"
	sortfile	= os.path.dirname(os.path.realpath(__file__)) + "/sort_file.pl"
	workdir  	= os.path.dirname(os.path.realpath(__file__))
	# example usage for help
	example  = "\tTo create a new project: specify (1) project name, (2) reference genomes and (3) input sequences (either --alignment, --fastq, or --mrsfast-best-search)\n"
	example += "\n\t--starting from a sam or bam file\n"
	example += "\t$ ./pamir.py -p my_project -r ref.fa --files alignment=my.sam\n"
	example += "\n\t--starting from a gzipped fastq file\n"
	example += "\t$ ./pamir.py -p my_project -r ref.fa --files fastq=my.fastq.gz\n"
	example += "\n\t--starting from a remapping result\n"
	example += "\t$ ./pamir.py -p my_project -r ref.fa --files mrsfast-best-search=my.best.sam\n"
	example += "\n\t--starting with a masked file\n"
	example += "\t$ ./pamir.py -p my_project -r ref.fa -m my-mask.txt  --files alignment=my.sam\n"	
	example += "\n\n\tTo resume a project, just type project folder and pamir.py will automatically resume from the previous stages:\n"
	example += "\t$ ./pamir.py -p my_project\n"
	example += "\t$ ./pamir.py -p /home/this/is/my/folder/project\n"


#############################################################################################
# Default values will be set up later in check_proj_preq
def command_line_process():
	parser = argparse.ArgumentParser(
		description='Pamir: Insertion Discovery in Whole Genome Sequencing Data',
		usage = pipeline.example,
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument('--project','-p',
		required=True,
		metavar='project',
		help='The name of the project. Pamir creates the folder if it does not exist'
	)
	parser.add_argument('--dest','-d',
		metavar='destination',
		help='Directory that will be used for analysis. (default: ./)'
	)
	parser.add_argument('--reference','-r',
		metavar='reference',
		help='The path to the reference genome that should be used for analysis'
	)
	parser.add_argument('--assembler','-a',
		metavar='assembler',
		help='The type of the assembler for assembling orphans: sga/minia. (default:sga)'
	)
	parser.add_argument('--max-contig','-c',
		type=int,
		metavar='max_contig',
		help="Maximum size of assembled contigs. (default: 400)",
	)
	parser.add_argument('--max-error','-e',
		type=int,
		metavar='max_error',
		help="Maximum number of mismatches within the structural variant. (default: 0)",
	)
	parser.add_argument('--min-support','-s',
		type=int,
		metavar='min_support',
		help="Support threshold (number of reads per base) for the predicted structural variants. (default: 5)",
	)
	parser.add_argument('--contig-identity','-i',
		type=int,
		metavar='contig_identity',
		help='Identity threshold (1..100) for an assembled contig to be considered for further analysis. (default: 85)',
	)
	parser.add_argument('--final-identity','-g',
		type=int,
		metavar='final_identity',
		help="Identity threshold (1..100) for the final alignment allowing for structural variants. (default: 95)",
	)	
	parser.add_argument('--mask-file', '-m',
		help='The coordinates provided in this file will be masked from the reference genome.'
	)
	parser.add_argument('--invert-masker',
		action='store_true',
		help='The provided coordinates will be masked in the reference genome and ignored in mapping step. (dafault: False)',
	)
	parser.add_argument('--mrsfast-index-ws',
		metavar='window_size',
		help='Window size used by mrsFAST-Ultra for indexing the reference genome. (default: 12)',
	)
	parser.add_argument('--mrsfast-errors',
		metavar='mrsfast_errors',
		help='Number of the errors allowed by mrsFAST-Ultra for mapping (default: -1)',
	)
	parser.add_argument('--mrsfast-threads',
		metavar='mrsfast_threads',
		help='Number of the threads used by mrsFAST-Ultra for mapping (default: 1)',
	)
	parser.add_argument('--mrsfast-n',
		metavar='mrsfast_n',
		help='Maximum number of mapping loci of anchor of an OEA. Anchor with higher mapping location will be ignored. 0 for considering all mapping locations (default: 1)',
	)
	parser.add_argument('--mrsfast-min',
		metavar='mrsfast_min',
		help='Minimum insert size for mrsFAST-Ultra paired-end best mapping (default: -1)',
	)
	parser.add_argument('--mrsfast-max',
		metavar='mrsfast_max',
		help='Maximum insert size for mrsFAST-Ultra paired-end best mapping (default: -1)',
	)
	parser.add_argument('--resume',
		nargs='?',
		const="pamir",
		help='Restart pipeline from the stage that has not been completed yet.',
	)
	parser.add_argument('--cluster',
		metavar='cluster',
		help='For input range x-y, report reads of cluster from x to y-1 and exit. Supported only in resume mode.',
	)
	parser.add_argument('--resume-force', '-f',
		action='store_true',
		help='Ignore existing files and restart the pipeline from stage specified in --resume.',
	)
	parser.add_argument('--num-worker',
		type=int,
		help='Number of independent prediction jobs to be created. (default: 1)',
	)
	parser.add_argument('--range',
		help='Intervals of OEA clusters in partition to be analyzed.',
	)
	parser.add_argument('--worker-id',
		help='Specific worker ID that the user wants to run prediction in the last stage. (default:-1, will run on all workers)',
	)
	parser.add_argument('--mode',
		metavar='engine_mode',
		help='Type for running indepedent jobs: normal, SGE, or PBS. (default: normal).',
		default='normal'
	)
	parser.add_argument('--job-max-time',
		help='Max job running time in PBS file. Read documents before changing its value! (default: 6 hour.)',
		default='06:00:00'
	)
	parser.add_argument('--job-max-memory',
		help='Max job memory in PBS file. Read documents before changing its value! (default: 16 GB.)',
		default='16G'
	)
	parser.add_argument('--files',
		metavar='files',
		nargs='+',
		type=str,
		default=[]
	)

	return parser

#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def log( msg, output=True ):
	if output:
		print "{0}".format(msg),
		sys.stdout.flush()
	
#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def logln( msg, output=True ):
	if output:
		print "{0}".format(msg)
		sys.stdout.flush()

#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def logOK():
	logln(bcolors.OKGREEN+"OK"+bcolors.ENDC)

#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def logWAR():
	logln(bcolors.WARNING+"SKIPPING"+bcolors.ENDC)

#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def logFAIL():
	logln(bcolors.FAIL+"FAILED"+bcolors.ENDC)

#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def shell(msg, run_command, command, control_file='', success_file='', success_msg='', shell_log=True):
	
	if shell_log:
		log("{0}...".format(msg))
	if run_command and command =='':
		print ""
	if (not run_command):
		logWAR()
		return 
	if not shell_log or control_file == '':
		p = subprocess.Popen(command, shell=True)
	else:
		f = open(control_file, 'w')
		f.write("CMD:" + command + "\n")
		f.flush() # so commands will show before actual progress
		p = subprocess.Popen(command, shell=True, stdout=f, stderr=f)

	ret = p.wait()
	if shell_log and control_file != '':
		f.close()
	if ret != 0:
		logFAIL()
		raise subprocess.CalledProcessError(ret, command)

	elif success_file != '':
		with open (success_file, 'w') as suc_file:
			suc_file.write(success_msg)
		logOK()
#############################################################################################
########## Clean stage file for resuming
def clean_state_worker( workdir, config):
	workdir = pipeline.workdir
	stage_dir = workdir + "/stage/"
	worker_stage_files = [f for f in os.listdir( stage_dir ) if (os.path.isfile( os.path.join( stage_dir , f)) and "30"==f[0:2])]
	for item in worker_stage_files:
		log(" + Removing " + item  + "...")
		os.remove( stage_dir  + item)
		logOK()
		
#############################################################################################
########## Clean stage file for resuming
def clean_state( mode_index, workdir, config ):
	flag_clean = 0 # 1 only when we delete files due to changes in project.config
	workdir = pipeline.workdir	
	#valid_state = [ '01.verify_sam', '02.mask', '03.mrsfast-index', '04.mrsfast.best', '05.remove_concordant', '06.mrsfast', '07.sorted', '08.oeaunm', '09.sniper_part', '10.orphan_assembly', '10.orphan_assembly_2','11.index_orphan_ref','16.oea2orphan','17.split1_oea2orphan','18.split2_oea2orphan','19.concat_all_oea2orphan','20.all_oea2orphan_recal','21.insert_orphans_into_cluster','22.updated_sniper_assembly','23.concatenate_vcf','24.index_log','25.sort_vcf','26.filter_duplicate_calls','27.remove_partials','28.generate_loclen','normal']
	valid_state = [ '01.verify_sam', '02.mask', '03.mrsfast-index', '04.mrsfast.best', '05.remove_concordant','normal']
	for i in range( mode_index, len(valid_state)):
		if ( 3==i and ""!=config.get("project","mrsfast-best-search")):
			continue
		if os.path.isfile(workdir + "/stage/" + valid_state[i] + ".finished" ):
		#try:
			if 0 == flag_clean:
				logln("Removing old files due to change in project.config")
				flag_clean = 1
			log(" + Removing " + valid_state[i] + ".finished...")
			os.remove( workdir + "/stage/" + valid_state[i] + ".finished" )
			logOK()
	# Clean all possible stages files of previous worker
	clean_state_worker( workdir, config)
	#except subprocess.CalledProcessError as e:
	#	print >>sys.stderr, "{0} failed with exit status {2} and message {1}".format(e.cmd, 'N/A', e.returncode)


#############################################################################################
###### Running commands for verify_sam 
def verify_sam(config ):
	msg           = "Sorting bam file"	
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/{1}".format(workdir, config.get("project","alignment"))
	output_file   = "{0}/{1}".format(workdir, config.get("project","fastq"))
	control_file  = "{0}/log/01.verify_sam_1.log".format(workdir);
	complete_file = "{0}/stage/01.verify_sam_1.finished".format(workdir);
	freeze_arg    = ""
	cmd			  = pipeline.samtools + " sort -n {0} {0}.sorted.bam".format(input_file) 
	run_cmd       = not (os.path.isfile(complete_file) )
	shell( msg, run_cmd , cmd, control_file, complete_file, freeze_arg)
	msg           = "Extracting FASTQ from Alignment file"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/{1}".format(workdir, config.get("project","alignment"))
	output_file   = "{0}/{1}".format(workdir, config.get("project","fastq"))
	control_file  = "{0}/log/01.verify_sam.log".format(workdir);
	complete_file = "{0}/stage/01.verify_sam.finished".format(workdir);
	freeze_arg    = ""
	cmd			  = pipeline.bedtools + "bamtofastq -i {0}.sorted.bam -fq {1}/input_1.fq -fq2 {1}/input_2.fq".format(input_file,output_file) 
#	cmd           = pipeline.sniper + ' verify_sam {0} {1}'.format( input_file, output_file )
	run_cmd       = not (os.path.isfile(complete_file) )
	if ( run_cmd ):
		clean_state( 1, workdir, config )
	shell( msg, run_cmd , cmd, control_file, complete_file, freeze_arg)
	
#############################################################################################
###### Running commands for mask 
def mask(config):
	msg			  = "Masking Reference Genome"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/{1}".format(workdir, config.get("project","reference"))
	mask_file     = "{0}/{1}".format(workdir, config.get("sniper","mask-file"))
	output_file   = "{0}/{1}.masked".format(workdir, config.get("project","reference"))
	mask_mode     = "maski" if ( "True" == config.get("sniper","invert-masker") ) else "mask"
	if ("maski" == mask_mode):
		msg += " to Keep Only Regions Specified in File"

	control_file  = "{0}/log/02.mask.log".format(workdir);
	complete_file = "{0}/stage/02.mask.finished".format(workdir);
	freeze_arg    = "invert" if "True" == config.get("sniper","invert-masker")  else "mask"
	cmd           = pipeline.sniper +' {0} {1} {2} {3} 0'.format( mask_mode, mask_file, input_file, output_file )
	run_cmd       = not ( os.path.isfile(complete_file) and freeze_arg in open(complete_file).read()) 

	if ( run_cmd ):
		clean_state( 2, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)

#############################################################################################
###### Running commands for each mode 
def index(config):
	msg           = "Indexing the masked genome using mrsFAST-Ultra"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/{1}.masked".format(workdir, config.get("project", "reference"))
	control_file  = "{0}/log/03.mrsfast-index.log".format(workdir);
	complete_file = "{0}/stage/03.mrsfast-index.finished".format(workdir);
	freeze_arg    = "ws={0}".format(config.get("mrsfast", "window_size") )
	cmd           = pipeline.mrsfast + ' --index  {0} --ws {1}'.format(input_file, config.get("mrsfast","window_size") )
	run_cmd       = not ( os.path.isfile(complete_file) and freeze_arg in open(complete_file).read()) 

	if ( run_cmd ):
		clean_state( 3, workdir, config)
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
	
	#return msg, run_cmd, cmd, control_file, complete_file, freeze_arg
#############################################################################################
###### Running commands for each mode 
def mrsfast_best_search(config):
	msg           = "Mapping non-concordant reads using mrsFAST-Ultra"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	index_file    = "{0}/{1}.masked".format(workdir, config.get("project", "reference"))
	input_file    = "{0}/{1}".format(workdir, config.get("project","fastq"))
	output_file   = "{0}/mrsfast.best.sam".format(workdir)
	threads       = config.get("mrsfast", "threads")
	control_file  = "{0}/log/04.mrsfast.best.log".format(workdir);
	complete_file = "{0}/stage/04.mrsfast.best.finished".format(workdir);
	freeze_arg    = "ws={0}.error={1}".format(config.get("mrsfast", "window_size"), config.get("mrsfast", "errors"))
	cmd           = pipeline.mrsfast +' --search {0} --threads {1} -n 0 --best --disable-sam-header --disable-nohits --pe --seq {2} --min {3} --max {4} --seqcomp -o {5}'.format(index_file, threads, input_file, config.get("mrsfast","min"), config.get("mrsfast","max"), output_file )
	if(config.get("mrsfast","errors")!="-1"):
		cmd+=" -e {0}".format(config.get("mrsfast","errors"))
	run_cmd       = not ( os.path.isfile(complete_file) and freeze_arg in open(complete_file).read()) 

	#if ( run_cmd ):
	#	clean_state( 4, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
	
	
#############################################################################################
###### Running commands for remove_concordant 
def remove_concordant(config,f):
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	msg           = "Extracting OEA and Orphan reads of {0}".format(f)
	input_file    = "{0}/mrsfast.best.sam".format(workdir)
	config.set("project","readlength", str(len((open(input_file).readline()).split("\t")[10])))
	with open ( workdir +"/project.config", "w") as configFile:
		config.write(configFile)
	output_file   = "{0}/".format(workdir)
	control_file  = "{0}/log/05.remove_concordant.log".format(workdir);
	complete_file = "{0}/stage/05.remove_concordant.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.sniper + ' remove_concordant {0} {1} 2 1 1'.format(  input_file, output_file )
	run_cmd       = not (os.path.isfile(complete_file) )
	if ( run_cmd ):
		clean_state( 5, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
#############################################################################################
###### Running commands for mrsfast_search
def mrsfast_search(config ):
	msg           = "Mapping OEA reads to get anchor locations using mrsFAST-Ultra"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	ref_file      = "{0}/{1}.masked".format(workdir, config.get("project", "reference"))
	input_file    = "{0}/oea.mapped.fq".format(workdir)
	output_file   = "{0}/mrsfast.sam".format(workdir)
	unmapped_file = "{0}/mrsfast.sam.unmap".format(workdir)
	threads       = config.get("mrsfast", "threads")
	cutoff        = config.get("mrsfast", "n")
	
	control_file  = "{0}/log/06.mrsfast.log".format(workdir);
	complete_file = "{0}/stage/06.mrsfast.finished".format(workdir);
	freeze_arg    = "ws={0}.error={1}.cutoff={2}".format(config.get("mrsfast", "window_size"), config.get("mrsfast", "errors"), cutoff)
	cmd           = pipeline.mrsfast +' --search  {0} --threads {1} --seq {2} -o {3} -u {4}'.format(ref_file, threads, input_file, output_file, unmapped_file )
	if(config.get("mrsfast","errors")!="-1"):
		cmd+=" -e {0}".format(config.get("mrsfast","errors"))
	if ("0"!= cutoff):
		cmd += " -n {0}".format(cutoff )
	run_cmd       = not ( os.path.isfile(complete_file) and freeze_arg in open(complete_file).read()) 

	if ( run_cmd ):
		clean_state( 6, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
#############################################################################################
###### Running commands for sorting all mapping output 
def sort(config ):
	msg           = "Sorting the mapped mates"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/mrsfast.sam".format( workdir )
	output_file   = "{0}/sorted.sam".format( workdir )
	control_file  = "{0}/log/07.sort.log".format( workdir );
	complete_file = "{0}/stage/07.sort.finished".format( workdir );
	freeze_arg    = ""
	cmd           = pipeline.sniper +' sort {0} {1}'.format(  input_file, output_file )
	run_cmd       = not (os.path.isfile(complete_file) )

	shell( msg, run_cmd , cmd, control_file, complete_file, freeze_arg)
#############################################################################################
###### Running commands for modifying oea.unmapped.fq  
def modify_oea_unmap(config ):
	msg           = "Modifying unmapped OEA read file format"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/oea.unmapped.fq".format(workdir )
	output_file   = "{0}/unmapped".format(workdir )
	control_file  = "{0}/log/08.modify_oea_unmap.log".format(workdir);
	complete_file = "{0}/stage/08.modify_oea_unmap.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.sniper + ' modify_oea_unmap {0} {1}'.format(  input_file, output_file )
	run_cmd       = not (os.path.isfile(complete_file) )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
#############################################################################################
###### Running commands for creating OEA clusters 
def sniper_part(config ):
	msg           = "Creating OEA clusters"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/sorted.sam".format(workdir)
	unmapped_file = "{0}/unmapped".format(workdir )
	output_file   = "{0}/sniper_part".format(workdir )
	control_file  = "{0}/log/09.sniper_part.log".format(workdir);
	complete_file = "{0}/stage/09.sniper_part.finished".format(workdir);
	freeze_arg    = "1000"
	cmd           = pipeline.sniper + ' partition {0} {1} {2} 1000'.format(  input_file, unmapped_file, output_file )
	run_cmd       = not (os.path.isfile(complete_file) )

	if ( run_cmd ):
		clean_state( 9, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
#############################################################################################
###### Assembling orphan reads into contigs 
def orphan_assembly(config):
	msg           = "Creating Orphan contigs"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	control_file  = "{0}/log/10.orphan_assembly.log".format(workdir);
	complete_file = "{0}/stage/10.orphan_assembly.finished".format(workdir);
	input_file    = "{0}/orphan.fq".format(workdir)
	freeze_arg    = ""
	msg = "Creating Orphan contigs"
	run_cmd       = not (os.path.isfile(complete_file) )
	#if ( run_cmd ):
	#	clean_state( 10, workdir, config )
	cmd=""
	if(config.get("project","assembler")=="velvet"):
		cmd 		  = pipeline.velveth + " {0}/velvetRes/ 31 -fastq -short {1}".format(workdir, input_file)
		shell(msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
		control_file  = "{0}/log/10.orphan_assembly_2.log".format(workdir);
		complete_file = "{0}/stage/10.orphan_assembly_2.finished".format(workdir);
		msg 		  = "Running velvetg"
		cmd 		  = pipeline.velvetg + " {0}/velvetRes/ -min_contig_lgth 400 -cov_cutoff 2".format(workdir)
		run_cmd       = not (os.path.isfile(complete_file) )
		shell(msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
		contigfile    = open("{0}/velvetRes/contigs.fa".format(workdir),"r").readlines()
		fout		  = open("{0}.contigs.fa".format(input_file),"w")
		a=0
		while a < len(contigfile):
			content =""
			if(contigfile[a][0]=='>'):
				conname = contigfile[a].strip().split("_")
				fout.write(">contig-"+conname[1]+"\n")
				a+=1
				while a < len(contigfile) and contigfile[a][0]!='>':
					content+=contigfile[a].strip()
					a+=1
				fout.write(content+"\n")
		fout.close()
	elif(config.get("project","assembler")=="sga"):
		cmd 		  = "python " + pipeline.sga + " " + input_file;
		shell(msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
	elif(config.get("project","assembler")=="minia"):
		cmd 		  = "{0} -in {1} -kmer-size 31 -abundance-min 3 -out {1}_tmp".format(pipeline.minia, input_file, input_file);
		shell(msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
		fin = open("{0}_tmp.contigs.fa".format(input_file),"r")
		fout = open("{0}.contigs.fa".format(input_file),"w")
		contigs = fin.readlines()
		i =0 
		while i < len(contigs):
			splitted=contigs[i].split("__")
			if int(splitted[2]) > 400:
				fout.write(">contig-"+splitted[0][1:]+ "\n" +contigs[i+1])
			i+=2
		fin.close()
		fout.close()
#############################################################################################
###### Preparing necessary file for orphan.contigs.fa to make it a single line ref: Easier to map.
def prepare_orphan_contig(config):
	msg			  = "Preparing orphan contigs for mapping"
	workdir		  = pipeline.workdir
	input_file    = "{0}/orphan.fq.contigs.fa".format(workdir)
	output_file   = "{0}/orphan.contigs.ref".format(workdir)
	output_file2  = "{0}/orphan.contigs.single.ref".format(workdir)
	coor_file     = "{0}/orphan.contigs.single.coor".format(workdir)
	f 			  = open (input_file, "r")
	fref 			  = open(output_file, "w")
	fsingleref 		  = open(output_file2, "w")
	fcoor			  = open(coor_file, "w")
	line=f.readline()
	fsingleref.write(">1\n")
	cloc=1;
	while(line!=""):
		cname= line.split(" ")[0]
		if cname.find('\n')==-1:
			cname=cname+'\n'
		fcoor.write(cname[1:len(cname)-1]+"\t"+str(cloc)+'\n')
		fref.write(cname)
		line=f.readline()
		cloc+=len(line)-1
		fref.write(line)
		fsingleref.write(line)
		line=f.readline()
	f.close()
	fref.close()
	fsingleref.close()
	fcoor.close()
	msg			  = "Indexing orphan.contigs.single.ref"
	cmd			  = pipeline.mrsfast + " --index {0}".format(output_file2)
	control_file  = "{0}/log/11.index_orphan_ref.log".format(workdir);
	complete_file = "{0}/stage/11.index_orphan_ref.finished".format(workdir);
	freeze_arg    = ""
	run_cmd       = not (os.path.isfile(complete_file) )
#	if ( run_cmd ):
#		clean_state( 11, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
#############################################################################################
###### Map orphan reads onto orphan contigs to find the support value of each contig (If they come from an assembler).
def orphan_to_orphan(config):
	msg			  = "Mapping orphan reads onto orphan contigs"
	workdir		  = pipeline.workdir
	orphan_ref    = "{0}/orphan.contigs.single.ref".format(workdir)
	orphan_fastq	  = "{0}/orphan.fq".format(workdir)
	output_file   = "{0}/orphan2orphan.sam".format(workdir)
	control_file  = "{0}/log/12.orphan2orphan.log".format(workdir);
	complete_file = "{0}/stage/12.orphan2orphan.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.mrsfast + ' --search {0} --seq {1} --crop {2} -o {3} -e 1 --disable-sam-header'.format(orphan_ref, orphan_fastq, config.get("project","readlength"),output_file)
	run_cmd       = not (os.path.isfile(complete_file) )
	if ( run_cmd ):
		clean_state( 12, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg )
	
	msg			  = "Recalibrating orphan_to_orphan.sam"
	input_file   = "{0}/orphan2orphan.sam".format(workdir)
	output_file   = "{0}/orphan2orphan.recal.sam".format(workdir)
	coor_file     = "{0}/orphan.contigs.single.coor".format(workdir)
	control_file  = "{0}/log/13.orphan2orphan.recal.log".format(workdir);
	complete_file = "{0}/stage/13.orphan2orphan.recal.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.recalibrate + " {0} {1} {2}".format(coor_file, input_file, output_file)
	freeze_arg=""
	run_cmd       = not (os.path.isfile(complete_file) )
	if ( run_cmd ):
		clean_state( 13, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg )
#############################################################################################
###### Map orphan reads onto orphan contigs to find the support value of each contig (If they come from an assembler).
def orphancontig_support(config):
	msg			  = "Sort orphan2orphan.recal.sam"
	workdir		  = pipeline.workdir
	input_file    = "{0}/orphan2orphan.recal.sam".format(workdir)
	output_file   = "{0}/orphan2orphan.recal.sorted.sam".format(workdir)
	control_file  = "{0}/log/14.sort_orphan_contig_sam.log".format(workdir);
	complete_file = "{0}/stage/14.sort_orphan_contig_sam.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.sniper + ' sort {0} {1}'.format(input_file, output_file)
	run_cmd       = not (os.path.isfile(complete_file) )
	if ( run_cmd ):
		clean_state( 14, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg )
	
	msg			  = "Printing contig-N_reads.fq for each orphan contig"
	input_file    = "{0}/orphan2orphan.recal.sorted.sam".format(workdir)
	control_file  = "{0}/log/15.orphan_contig_support.log".format(workdir);
	complete_file = "{0}/stage/15.orphan_contig_support.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.ext_sup + " {0} {1}".format(input_file, pipeline.workdir)
	freeze_arg=""
	run_cmd       = not (os.path.isfile(complete_file) )
	if ( run_cmd ):
		clean_state( 15, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg )
#############################################################################################
###### Map oea reads onto orphan contigs.
def oea_to_orphan(config):
	msg			  = "Mapping oea reads onto orphan contigs"
	workdir		  = pipeline.workdir
	orphan_ref    = "{0}/orphan.contigs.single.ref".format(workdir)
	seq_fastq	  = "{0}/oea.unmapped.fq".format(workdir)
	output_file   = "{0}/oea2orphan.sam".format(workdir)
	control_file  = "{0}/log/16.oea2orphan.log".format(workdir);
	complete_file = "{0}/stage/16.oea2orphan.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.mrsfast + ' --search {0} --seq {1} -o {3} -e 0 --disable-sam-header'.format(orphan_ref, seq_fastq, config.get("project","readlength"),output_file)
	run_cmd       = not (os.path.isfile(complete_file) )
#	if ( run_cmd ):
#		clean_state( 16, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg )
	
################################################################################################
###### Map nohit oea reads as split reads onto orphan contigs.
def oea_to_orphan_split(config):
	msg			  = "Mapping first readlength/2 of OEA reads as split onto orphan contigs"
	workdir		  = pipeline.workdir
	orphan_ref    = "{0}/orphan.contigs.single.ref".format(workdir)
	input_fastq	  = "{0}/oea2orphan.sam.nohit".format(workdir)
	output_file   = "{0}/split1_oea2orphan.sam".format(workdir)
	control_file  = "{0}/log/17.split1_oea2orphan.log".format(workdir);
	complete_file = "{0}/stage/17.split1_oea2orphan.finished".format(workdir);
	freeze_arg    = ""
	half 		  = (int(config.get("project","readlength"))/2)
	cmd           = pipeline.mrsfast + ' --search {0} --seq {1} -o {2} --crop {3} -e 0 --disable-sam-header'.format(orphan_ref, input_fastq, output_file, str(half))
	run_cmd       = not (os.path.isfile(complete_file) )
#	if ( run_cmd ):
#		clean_state( 17, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg )


	msg			  = "Mapping last readlength/2 of OEA reads as split onto orphan contigs"
	workdir		  = pipeline.workdir
	orphan_ref    = "{0}/orphan.contigs.single.ref".format(workdir)
	input_fastq	  = "{0}/oea2orphan.sam.nohit".format(workdir)
	output_file   = "{0}/split2_oea2orphan.sam".format(workdir)
	control_file  = "{0}/log/18.split2_oea2orphan.log".format(workdir);
	complete_file = "{0}/stage/18.split2_oea2orphan.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.mrsfast + ' --search {0} --seq {1} -o {2} --tail-crop {3} -e 0 --disable-sam-header'.format(orphan_ref, input_fastq, output_file, str(half))
	run_cmd       = not (os.path.isfile(complete_file) )
#	if ( run_cmd ):
#		clean_state( 18, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg )
	
	msg			   = "Concatenating all oea2orphan"
	input_file1    = "{0}/split1_oea2orphan.sam".format(workdir)
	input_file2    = "{0}/split2_oea2orphan.sam".format(workdir)
	input_file3    = "{0}/oea2orphan.sam".format(workdir)
	output_file    = "{0}/all_oea2orphan.sam".format(workdir)
	control_file   = "{0}/log/19.concat_all_oea2orphan.log".format(workdir);
	complete_file  = "{0}/stage/19.concat_all_oea2orphan.finished".format(workdir);
	freeze_arg     = ""
	cmd            = "cat {0} {1} {2} > {3}".format(input_file1, input_file2, input_file3, output_file)
	freeze_arg=""
	run_cmd        = not (os.path.isfile(complete_file) )
	#if ( run_cmd ):
	#	clean_state( 19, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg )
#############################################################################################
###### Recalibrate all oea2orhan.sam.
def recalibrate_all_oea_to_orphan(config):
	msg			  = "Recalibrating the all_oea_to_orphan.sam"
	workdir		  = pipeline.workdir
	input_file    = "{0}/all_oea2orphan.sam".format(workdir)
	output_file   = "{0}/all_oea2orphan.recal.sam".format(workdir)
	coor_file     = "{0}/orphan.contigs.single.coor".format(workdir)
	control_file  = "{0}/log/20.all_oea2orphan_recal.log".format(workdir);
	complete_file = "{0}/stage/20.all_oea2orphan_recal.finished".format(workdir);
	freeze_arg    = ""
	cmd           = "{0} {1} {2} {3}".format(pipeline.recalibrate, coor_file, input_file, output_file)
	freeze_arg=""
	run_cmd       = not (os.path.isfile(complete_file) )
#	if ( run_cmd ):
#		clean_state( 20, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg )
#############################################################################################
###### Insert orphan contigs into OEA clusters.
def orphans_into_oeacluster(config):
	msg			  = "Inserting relevant orphan contigs into oea cluster"
	workdir		  = pipeline.workdir
	orphan_ref    = "{0}/orphan.contigs.ref".format(workdir)
	orphan_to_orphan	  = "{0}/orphan2orphan.recal.sam".format(workdir)
	oea_to_orphan   = "{0}/all_oea2orphan.recal.sam".format(workdir)
	partition_file   = "{0}/sniper_part".format(workdir)
	upartition_file   = "{0}/sniper_part_updated".format(workdir)
	control_file  = "{0}/log/21.insert_orphans_into_cluster.log".format(workdir)
	complete_file = "{0}/stage/21.insert_orphans_into_cluster.finished".format(workdir)
	clusterNumFile = "{0}/clusterNum".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.pprocessor + ' {0} {1} {2} {3} {4} > {5}'.format(orphan_ref, orphan_to_orphan, oea_to_orphan, partition_file, upartition_file, clusterNumFile)
	run_cmd       = not (os.path.isfile(complete_file) )
	if ( run_cmd ):
		clean_state( 21, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg )
#############################################################################################
###### Print VCF header 
def print_header(config):
	msg           = "Printing VCF header"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	header_file    = "{0}/vcf_header".format(workdir)
	ref_file	  = "{0}/{1}".format(workdir,config.get("project","reference"))
	control_file  = "{0}/log/22.printheader.log".format(workdir)
	complete_file = "{0}/stage/22.printheader.finished".format(workdir)
	cmd           = pipeline.sniper + ' header {0} {1}'.format( header_file, ref_file)
	run_cmd       = not (os.path.isfile(complete_file))
	if ( run_cmd ):
		clean_state( 22, workdir, config )
	freeze_arg    = ""
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
######################################################################################
###### Assemble updated cluster 
def updated_sniper_part(config ):
	msg           = "Assembling updated cluster"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/sniper_part_updated".format(workdir)
	ref_file	  = "{0}/{1}".format(workdir,config.get("project","reference"))
	control_file  = "{0}/log/23.updated_sniper_assembly.log".format(workdir)
	complete_file = "{0}/stage/23.updated_sniper_assembly.finished".format(workdir)
	clusterNumFile= "{0}/clusterNum".format(workdir)
	clusterNum=int(open(clusterNumFile).read())
	perJob=clusterNum/int(config.get("project","num-worker"))
	jobFileName= "{0}/run_sniper_assemble.sh".format(workdir)
	f     = open("{0}".format(jobFileName),"w")
	f.write("#!/bin/bash\n")
	i=0
	workNum=1
	while i<clusterNum:
		output_file   = "{0}/sniper_part_updated.vcfx{1}".format(workdir,workNum)
		output_log	  = "{0}/sniper_part_updated.logx{1}".format(workdir,workNum)
		cmd           = pipeline.sniper + ' assemble {0} {1} {2}-{3} {4} 30000 {5} {7} > {6}'.format( input_file, ref_file, str(i),str(i+perJob), output_file, str(config.get("project","readlength")), output_log, workdir)
		f.write(cmd+'\n')
		i+=perJob
		workNum+=1
	f.close()
	freeze_arg    = ""
	cmd           =  'cat {0} | xargs -I CMD --max-procs={1} bash -c CMD'.format( jobFileName, config.get("project","num-worker"))
	run_cmd       = not (os.path.isfile(complete_file) )
	if ( run_cmd ):
		clean_state( 23, workdir, config )
	
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)

### concatenate outputs
	i=1
	msg = "Concatenating all vcf and files"
	cmd="cat "
	cmd2="cat "
	cmd3="cat "
	while i< workNum:
		cmd+="{0}/sniper_part_updated.vcfx{1} ".format(workdir,i)
		cmd2+="{0}/sniper_part_updated.logx{1} ".format(workdir,i)
		i+=1
	cmd+="> {0}/sniper_part_updated.vcf".format(workdir)
	cmd2+="> {0}/sniper_part_updated.log".format(workdir)
	cmdall=cmd+";"+cmd2
	freeze_arg=""
	control_file  = "{0}/log/24.concatenate_vcf.log".format(workdir)
	complete_file = "{0}/stage/24.concatenate_vcf.finished".format(workdir)
	run_cmd       = not (os.path.isfile(complete_file) )
	shell( msg, run_cmd, cmdall, control_file, complete_file, freeze_arg)

### index sniper_part_updated.log

	msg = "Indexing log file"
	cmd = pipeline.sniper + " index_log {0}/sniper_part_updated.log".format(workdir)
	freeze_arg=""
	control_file  = "{0}/log/25.index_log.log".format(workdir)
	complete_file = "{0}/stage/25.index_log.finished".format(workdir)
	run_cmd       = not (os.path.isfile(complete_file) )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
######################################################################################
#### sort the vcf and remove the duplications generate interleaved and paired fastq files from orphans and unmapped oeas.
def dupremoval_cleaning(config):	
	workdir		  = pipeline.workdir
	freeze_arg=""
	control_file  = "{0}/log/27.sort_and_filter_duplicate_calls.log".format(workdir)
	complete_file = "{0}/stage/27.sort_and_filter_duplicate_calls.finished".format(workdir)
	run_cmd		  = not (os.path.isfile(complete_file))
	cmd			  = pipeline.sortvcf + " {0}/sniper_part_updated.vcf {0}/sniper_part_updated.vcf.sorted_wodups 1".format(workdir)
	msg="Sorting VCF file and eliminating duplicated insertions"
	shell(msg,run_cmd,cmd,control_file,complete_file,freeze_arg)
	msg = "You can check output file now: sniper_part_updated.vcf.sorted_wodups"
	shell(msg,True,"")
	control_file  = "{0}/log/28.remove_partials.log".format(workdir)
	complete_file = "{0}/stage/28.remove_partials.finished".format(workdir)
	run_cmd       = not (os.path.isfile(complete_file))
	cmd="rm {0}/sniper_part_updated.vcfx* {0}/sniper_part_updated.logx*".format(workdir)
	msg="Deleting partial outputs"
	shell(msg,run_cmd,cmd,control_file,complete_file,freeze_arg)
	control_file  = "{0}/log/29.generate_loclen.log".format(workdir)
	complete_file = "{0}/stage/29.generate_loclen.finished".format(workdir)
	run_cmd       = not (os.path.isfile(complete_file))
	cmd="cut -f2,3 {0}/sniper_part_updated.vcf.sorted_wodups > {0}/sniper_part_updated.vcf.sorted_wodups_loclen".format(workdir)
	msg="Generating _loclen"
	shell(msg,run_cmd,cmd,control_file,complete_file,freeze_arg)
######################################################################################
#### filtering and genotyping.
def post_processing(config):	
	workdir		  = pipeline.workdir
	TLEN = 1000
	freeze_arg=""
	control_file  = "{0}/log/30.filtering.log".format(workdir)
	complete_file = "{0}/stage/30.filtering.finished".format(workdir)
	run_cmd		  = not (os.path.isfile(complete_file))
	cmd			  = pipeline.filtering + " {0}/sniper_part_updated.vcf.sorted_wodups {0}/{1}.masked {2} {3} {4} {0} {5}".format(workdir, config.get("project","reference"), config.get("project","readlength"), config.get("mrsfast","min"), config.get("mrsfast","max"), str(TLEN))
	msg="Filtering insertion candidates"
	shell(msg,run_cmd,cmd,control_file,complete_file,freeze_arg)
	###### Prepare input file for setcover
	control_file  = "{0}/log/31.generate_set_cover_input.log".format(workdir)
	complete_file = "{0}/stage/31.generate_set_cover_input.finished".format(workdir)
	run_cmd		  = not (os.path.isfile(complete_file))
	cmd			  = pipeline.gensetcov + " {0}/sniper_part_updated.vcf.sorted_wodups_filtered_forSETCOVER.sorted {0}/filtering/seq.mrsfast.recal.sam.sorted {0}/forSETCOVER".format(workdir)
	msg="Preparing input file for setcover"
	shell(msg,run_cmd,cmd,control_file,complete_file,freeze_arg)
	###### Run setcover
	control_file  = "{0}/log/32.setcover.log".format(workdir)
	complete_file = "{0}/stage/32.setcover.finished".format(workdir)
	run_cmd		  = not (os.path.isfile(complete_file))
	cmd			  = pipeline.smoother + " {0}/forSETCOVER > {0}/fromSETCOVER".format(workdir)
	msg="Running setcover"
	shell(msg,run_cmd,cmd,control_file,complete_file,freeze_arg)
	###### Eliminate the ones removed by setcover
	control_file  = "{0}/log/33.filter_by_setcover.log".format(workdir)
	complete_file = "{0}/stage/33.filter_by_setcover.finished".format(workdir)
	run_cmd		  = not (os.path.isfile(complete_file))
	cmd			  = pipeline.filterbysetcover + " {0}/fromSETCOVER {0}/sniper_part_updated.vcf.sorted_wodups_filtered {0}/sniper_part_updated.vcf.sorted_wodups_filtered_aftersetcov".format(workdir)
	msg="Filter removed calls by setcover"
	shell(msg,run_cmd,cmd,control_file,complete_file,freeze_arg)
	###### Grep PASS calls
	control_file  = "{0}/log/34.grepPASS.log".format(workdir)
	complete_file = "{0}/stage/34.grepPASS.finished".format(workdir)
	run_cmd		  = not (os.path.isfile(complete_file))
	cmd			  = "grep PASS {0}/sniper_part_updated.vcf.sorted_wodups_filtered_aftersetcov > {0}/sniper_part_updated.vcf.sorted_wodups_filtered_aftersetcov_PASS".format(workdir)
	msg="Grep PASS calls"
	shell(msg,run_cmd,cmd,control_file,complete_file,freeze_arg)
	###### Sort setcover output
	control_file  = "{0}/log/35.sortSETCover.log".format(workdir)
	complete_file = "{0}/stage/35.sortSETCover.finished".format(workdir)
	run_cmd		  = not (os.path.isfile(complete_file))
	cmd			  = "perl " + pipeline.sortfile + " {0}/sniper_part_updated.vcf.sorted_wodups_filtered_aftersetcov_PASS".format(workdir)
	msg="Sort setcover output"
	shell(msg,run_cmd,cmd,control_file,complete_file,freeze_arg)
#############################################################################################
###### Running commands for extracting clusters 
def output_cluster(config, c_range ):
	msg           = "Extracting Clusters for range {0}".format(c_range)
	workdir		  = pipeline.workdir
	control_file  = "{0}/log/output_cluster.log".format(workdir);
	cmd           = pipeline.sniper + " get_cluster {0}/sniper_part {1}".format(workdir, c_range)
	shell( msg, True, cmd, control_file, '', '')
#############################################################################################
def remove_concordant_for_each_bestsam(config):
	workdir = pipeline.workdir
	inputfiles = os.listdir(workdir+'/bestsam')
	if not os.path.isfile("{0}/stage/05.all_remove_concordant.finished".format(workdir)):
		cmd = "cat "
		cmd2 = "cat "
		cmd3 = "cat "
		cmd4 = "cat "
		a=1
		for f in inputfiles:
			if os.path.isfile("{0}/stage/05.remove_concordant.finished".format(workdir)):
				os.system("rm {0}/stage/05.remove_concordant.finished".format(workdir))
			if os.path.isfile("{0}/mrsfast.best.sam".format(workdir)):
				os.system("rm {0}/mrsfast.best.sam".format(workdir))
			if os.path.isfile(workdir+"/bestsam/"+f):
				symlink_name(workdir+"/bestsam/"+f, workdir, "mrsfast.best.sam")
			config.set("project", "mrsfast-best-search", os.path.abspath(f))
			remove_concordant(config,f)
			os.system("mv {0}/orphan.fq {0}/{1}.orphan.fq".format(workdir,a))
			os.system("mv {0}/oea.mapped.fq {0}/{1}.oea.mapped.fq".format(workdir,a))
			os.system("mv {0}/oea.unmapped.fq {0}/{1}.oea.unmapped.fq".format(workdir,a))
			os.system("mv {0}/all_interleaved.fastq {0}/{1}.all_interleaved.fastq".format(workdir,a))
			cmd+="{0}/{1}.orphan.fq ".format(workdir,a)
			cmd2+="{0}/{1}.oea.mapped.fq ".format(workdir,a)
			cmd3+="{0}/{1}.oea.unmapped.fq ".format(workdir,a)
			cmd4+="{0}/{1}.all_interleaved.fastq ".format(workdir,a)
			a+=1
		cmd+= " > {0}/orphan.fq".format(workdir)
		cmd2+= " > {0}/oea.mapped.fq".format(workdir)
		cmd3+= " > {0}/oea.unmapped.fq".format(workdir)
		cmd4+= " > {0}/all_interleaved.fastq".format(workdir)
		os.system(cmd)
		os.system(cmd2)
		os.system(cmd3)
		os.system(cmd4)
		os.system("touch {0}/stage/05.all_remove_concordant.finished".format(workdir))

#############################################################################################
###### Running commands for each mode 
def run_command(config, force=False):
	verify_sam(config)
	mask(config)
	index(config)
	mrsfast_best_search(config)
	remove_concordant_for_each_bestsam(config)
	mrsfast_search(config)
	sort(config)
	modify_oea_unmap(config)
	sniper_part(config)
	orphan_assembly(config)
	prepare_orphan_contig(config)
	oea_to_orphan(config)
	oea_to_orphan_split(config)
	recalibrate_all_oea_to_orphan(config)
	orphans_into_oeacluster(config)
	print_header(config)
	updated_sniper_part(config)
	exit(1)
	dupremoval_cleaning(config)
	post_processing(config)
	exit(0)

#############################################################################################
def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as e:
		if e.errno == errno.EEXIST and os.path.isdir(path):
			print "[ERROR] The project folder exists. Please run in resume mode or delete the project folder to re-run from scratch"
			exit(1);
#############################################################################################
######### link the absolute path of src in dest with identical filename
def symlink(src, dest):
	if src == '':
		return
	if not os.path.isfile(src) and not os.path.isdir(src):
		print "[ERROR] Input file {0} does not exist".format(src)
		exit(1)
	basename = os.path.basename(src)
	dest=os.path.abspath(dest)+'/'+basename
	#abs_src= os.path.abspath(src) + "/" + basename
	#os.symlink(src, dest)
	os.symlink(os.path.abspath(src), dest)

#############################################################################################
######### link the absolute path of src in dest with new filename
def symlink_name(src, dest, filename ):
	if src == '':
		return
	if not os.path.isfile(src) and not os.path.isdir(src):
		print "[ERROR] Input file {0} does not exist".format(src)
		exit(1)
	basename = os.path.basename(filename)
	dest=os.path.abspath(dest)+'/'+basename
	#abs_src= os.path.abspath(src) + "/" + basename
	#os.symlink(src, dest)
	os.symlink(os.path.abspath(src), dest)
#############################################################################################
##########	Make sure the project is successfully built
def is_exec(f):
	return os.path.isfile(f) and os.access(f, os.X_OK)
#############################################################################################
def check_binary_preq():
	execs = ['mrsfast', 'sniper']
	log( "Checking binary pre-requisites... ")
	for exe in execs:
		exe = os.path.dirname(os.path.realpath(__file__)) + "/" + exe
		if not is_exec(exe):
			print "[ERROR] File {0} cannot be executed. Please use 'make' to build the required binaries.".format(exe)
			logFAIL()
			logln ("File {0} cannot be executed. Please use 'make' to build the required binaries.".format(exe) )
			exit(1)
		
	logOK()

#############################################################################################
def resume_state_help():
	print "\nPamir supports the following resume states:"
	print "\tverify_sam: extract fastq file from given alignment"
	print "\tmask: masking reference genome"
	print "\tmrsfast-index: indexing reference genome for further mapping"
	print "\tmrsfast-best-search: mapping candidate reads"
	print "\tremove_concordant: extract oea reads"
	print "\tmrsfast-search: mapping oea reads to collect anchor locations"
	print "\tsort: sorting reads according to anchor locations"
	print "\tmodify_oea_unmap: extract unmapped mates for predicting microSVs"
	print "\tsniper_part: cluster unmapped reads according to anchor locations"
	print "\tnum-worker: generate jobs for parallel processing"
	print "\tsniper: output micorSVs"
	print "\nNOTE\tIf you want to automatically resume a killed job, just type --resume"

#############################################################################################
# Checking if necessary files are provided for a NEW project.
# reference sequence
# at least sam or fastq needs to exist
# masking file: if not provided, created a empty file
def check_input_preq( config ):
	workdir = pipeline.workdir

	if (None == config.get("project", "reference") ) or (not os.path.isfile( config.get("project", "reference") ) ):
		logFAIL()
		logln("Reference Genome {0} does not exist. Please provide valid name or path.".format( config.get("project", "reference") ))
		exit(1)
	
	if ("" !=  config.get("sniper", "mask-file") and  not os.path.isfile( config.get("sniper", "mask-file") )):
		logFAIL()
		logln("Invalid mask file {0}. Please provide a correct path.".format( os.path.abspath( config.get("sniper", "mask-file")) ))
		exit(1)

	inp = False # detect any valid input file or not
	if "" != config.get("project", "mrsfast-best-search"):
		inp = True
		if not (os.path.isfile(config.get("project", "mrsfast-best-search")) or os.path.isdir(config.get("project","mrsfast-best-search")) or config.get("project","mrsfast-best-search").find(",")==-1) :
			logFAIL()
			logln("mrsFAST remapping file {0} does not exist. Please provide a correct path.".format( config.get("project", "mrsfast-best-search")))
			exit(1)
	if "" != config.get("project", "fastq"):
		if (inp) :
			logFAIL()
			logln("Please provide either alignment file, fastq file, or best-search file.")
			exit(1)
		else:
			inp = True
		if not os.path.isfile(config.get("project", "fastq")):
			logFAIL()
			logln("Fastq file {0} does not exist. Please provide a correct path.".format( config.get("project", "fastq")))
			exit(1)
	
	if "" != config.get("project", "alignment"):
		if (inp) :
			logFAIL()
			logln("Please provide either alignment file, fastq file, or best-search file.")
			exit(1)
		else:
			inp = True
		if not os.path.isfile(config.get("project", "alignment")):
			logFAIL()
			logln("Alignment file {0} does not exist. Please provide a correct path.".format( config.get("project", "alignment")))
			exit(1)
	
	if not inp:
		logFAIL()
		logln("Please provide either alignment file, fastq file, or best-sesarch file.")
		exit(1)

#############################################################################################
def get_input_file(config, args_files):
	for i in args_files:
		i = i.split('=')
		if i[0]=="alignment": 
			config.set("project","alignment",i[1])
			#start_from_fastq = 0 # start from alignment
		elif i[0]=='fastq':
			config.set("project","fastq",i[1])
			#start_from_fastq = 1 # start from fastq
		elif i[0]=='mask':
			config.set("project",'mask',i[1])
			#start_from_fastq = 1 # start from fastq
		elif i[0]=='mrsfast-index':
			config.set("project",'mrsfast-index',i[1])
			#start_from_fastq = 1 # start from fastq
		elif i[0]=='mrsfast-best-search':
			config.set("project",'mrsfast-best-search',i[1])
		#elif i[0]=='oea':
		#	config.set("project",'oea',i[1])
		#elif i[0]=='mrsfast-search':
		#	config.set("project",'mrsfast-search',i[1])
		#elif i[0]=='sort':
		#	config.set("project",'sort',i[1])
		#elif i[0]=='oeanum':
		#	config.set("project",'oeanum',i[1])
	return config

#############################################################################################
########## Initialze mrsfast parameters for before creating project folder
def initialize_config_mrsfast( config, args):
	config.add_section("mrsfast")
	config.set("mrsfast", "window_size", str( args.mrsfast_index_ws ) if args.mrsfast_index_ws != None else  "12")
	config.set("mrsfast", "threads", str( args.mrsfast_threads ) if args.mrsfast_threads != None else  "1")
	config.set("mrsfast", "errors", str( args.mrsfast_errors ) if args.mrsfast_errors != None else  "-1")
	config.set("mrsfast", "min", str( args.mrsfast_min ) if args.mrsfast_min != None else  "-1")
	config.set("mrsfast", "max", str( args.mrsfast_max ) if args.mrsfast_max != None else  "-1")
	config.set("mrsfast", "n", str( args.mrsfast_n ) if args.mrsfast_n != None else  "50")
	return config

#############################################################################################
########## Initialze sniper parameters for before creating project folder
def initialize_config_sniper( config, args):
	workdir = pipeline.workdir
	config.add_section("sniper")
	config.set("sniper", "max-contig", str( args.max_contig ) if args.max_contig != None else  "25000")
	config.set("sniper", "max-error", str( args.max_error ) if args.max_error != None else "1")
	config.set("sniper", "mask-file", args.mask_file if args.mask_file != None else "" )
	config.set("sniper","engine-mode",args.mode if args.mode != None else "normal")
	config.set("sniper","invert-masker", str(args.invert_masker) if args.invert_masker !=None else "False")
	config.set("sniper","range", str( args.range ) if args.range !=None else "-1" )
	config.set("sniper","worker-id", str(args.worker_id) if args.worker_id != None else "-1")
	config.set("sniper","job-time", args.job_max_time if args.job_max_time != None else "06:00:00")
	config.set("sniper","job-memory",args.job_max_memory if args.job_max_memory != None else "16G")
	return config

#############################################################################################
def check_project_preq():
	args = command_line_process().parse_args()
	config = ConfigParser.ConfigParser()
	
	project_name = os.path.basename(os.path.normpath(args.project))
	pipeline.workdir = os.path.abspath(args.project)
	workdir = pipeline.workdir

	print "============================================="
	print "Project Name      : "+bcolors.OKGREEN+project_name+bcolors.ENDC
	print "Working Directory : "+bcolors.OKGREEN+workdir+bcolors.ENDC
	print "============================================="
	
	# Check if users want to resume the project
	if ( os.path.isdir( workdir)):
		log ("Checking the project pre-requisites... ")
		if ( None == args.resume):
			logFAIL()
			logln("Pamir can not overwrite an existing project. Please add --resume or change project name.")
			exit(1)
			
		if not os.path.isfile( workdir + "/project.config"):
			logFAIL()
			logln("NO config settings found. Please remove and re-create the project.")
			exit(1)
		logOK()
		
		log ("Loading the config file... ")
		config.read(workdir + '/project.config');
		# update range  and worker id for assemble stage in SGE and PBS
		config.set("sniper","range", str( args.range ) if args.range !=None else "-1" )
		config.set("sniper","worker-id", str( args.worker_id ) if args.worker_id !=None else "-1" )
		logOK()
		# Entering cluster extraction mode
		if ( args.cluster != None ):
			output_cluster( config, args.cluster )	
			logOK()
			exit(0)

	# Start a new project
	else:
		log("Creating a new project folder...")
		# set up main project parameter
		config.add_section("project")
		config.set("project", "name", project_name)
		config.set("project", "reference", args.reference)
		config.set("project", "num-worker", str(args.num_worker) if args.num_worker != None else "1" )
		config.set("project", "alignment",'')
		config.set("project", "fastq",'')
		config.set("project", "mrsfast-best-search",'')
		config.set("project", "assembler",str(args.assembler) if args.assembler!=None else "velvet")
		config = get_input_file( config, args.files)

		# Parameters for other parts in the pipeline
		initialize_config_mrsfast(config, args)
		initialize_config_sniper(config, args)
		
		#validating required files according to mode
		check_input_preq(config)
		
		# creating project folder
		mkdir_p(workdir)
		mkdir_p(workdir +'/jobs');
		mkdir_p(workdir +'/log');
		mkdir_p(workdir +'/pbs');
		mkdir_p(workdir +'/stage');

		# linking reference genome
		symlink(config.get("project", "reference"), workdir)
		config.set("project", "reference", os.path.basename(config.get("project", "reference")))

		if ("" != config.get("sniper", "mask-file")):
			symlink(config.get("sniper", "mask-file"),  workdir)
			config.set("sniper", "mask-file", os.path.basename(config.get("sniper", "mask-file")))

		# Exactly one non-empty path for input sources. set up files based on input source
		if ( "" != config.get("project", "alignment")):
			symlink(config.get("project", "alignment"), workdir)
			config.set("project", "alignment", os.path.basename(config.get("project", "alignment")) )
			config.set("project", "fastq",  project_name + ".fastq" )
			
		elif( "" != config.get("project", "fastq")):
			symlink(config.get("project", "fastq"), workdir)
			config.set("project", "fastq", os.path.basename(config.get("project", "fastq")))
		elif( "" != config.get("project", "mrsfast-best-search")):
			if(os.path.isdir(config.get("project","mrsfast-best-search"))):
				symlink_name(config.get("project","mrsfast-best-search"),workdir,'/bestsam')
			else:
				mkdir_p(workdir + '/bestsam');
				flist =  config.get("project","mrsfast-best-search").split(",")
				for f in flist:
					symlink_name(f,workdir+'/bestsam/',f)

		#	symlink_name(config.get("project", "mrsfast-best-search"), workdir, "mrsfast.best.sam")
		#	config.set("project", "mrsfast-best-search", os.path.basename(config.get("project", "mrsfast-best-search")))
			# Make sure mrsfast-best-search will not start when giving best mapping file
			with open("{0}/stage/04.mrsfast.best.finished".format(workdir), 'w') as complete_file:
				complete_file.write("ws={0}.error={1}\n".format(config.get("mrsfast", "window_size"), config.get("mrsfast", "errors")) )
		#else:
		#	loglng("[ERROR] Multiple input data sources. Check --files parameters before re-run")
		#	exit(1)

		# creating config in folder
		with open ( workdir +"/project.config", "w") as configFile:
			config.write(configFile)
		logOK()


	if ("" == config.get("sniper", "mask-file")):
		file = open ("{0}/mask.txt".format(workdir),'w')
		file.close()
		config.set("sniper","mask-file", "mask.txt")

	# set up stage files. Note that exactly one of three files are non-empty now
	if ("" == config.get("project", "alignment") ):
		file = open ("{0}/stage/01.verify_sam.finished".format(workdir), 'w')
		file.close()
		file = open ("{0}/stage/01.verify_sam_1.finished".format(workdir), 'w')

	return config
#############################################################################################
def main():
	config = check_project_preq()
	check_binary_preq()

	resume_state="sniper"

	try:
		if config.get("sniper", "range") == '-1':
			run_command(config)
		elif resume_state == "sniper":
			assemble(config)
		else:
			raise Exception('Invalid mode selected: ' + mode)
	except subprocess.CalledProcessError as e:
		print >>sys.stderr, "{0} failed with exit status {2} and message {1}".format(e.cmd, 'N/A', e.returncode)

#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

