#!/bin/bash
#set -e
echo "Running Pamir ($(git rev-parse  --short HEAD))"
SCRIPT_PATH="$(dirname `which $0`)/pamir"

#command -v conda && echo "Found conda, skipping checks :)" >&2 || { NO_CONDA=1; }

VERS_CHECKER="$SCRIPT_PATH/scripts/version_check.py"


ex(){
    command -v $1 1> /dev/null 2> /dev/null

    if [ "$?" -eq "1" ]
    then
        echo "$1... Failed, Not in PATH"
        exit 1
    fi
}

vx(){
    printf "$1($2)..."

    if [ "$3" -eq 0 ]
    then
        echo " OK, Version not checked."
        return
    fi
    vers_checker_str=$(cat environment.yaml | grep $1 | awk '{print $(NF-1)" "$(NF)}')
    

    python $VERS_CHECKER  $vers_checker_str $2


    
    if [ $? -eq "0" ]
    then
        echo " OK "
    else
        echo " Failed, requiring $vers_checker_str"
        exit 1
    fi


}

CFG_PATH=$(echo $@ | sed 's/=/ /g' |awk '{for(i=1;i <= NF;i++){if($i=="--configfile"){print $(i+1)}}}')

ex "snakemake" 
vx "snakemake" $(snakemake --version |head -n 1| sed 's/[A-Za-z\+]//g') 1


ex "samtools" 
vx "samtools" $(samtools  --version |head -n 1| awk '{print $2;}' | sed 's/[A-Za-z\+]//g') 1

ex "bedtools" 
vx "bedtools" $( bedtools  --version |head -n 1| awk '{print $2;}' | cut -d"-" -f1 | sed 's/[A-Za-z\+]//g') 1

ex "mrsfast"
vx "mrsfast" $( mrsfast   --version |head -n 1| awk '{print $2;}' | sed 's/[A-Za-z\+]//g') 1

ex "bwa"
vx "bwa" $(bwa 2>&1 | grep Version | awk '{print $2}' | awk -F"-" '{print $1;}') 1

ex "RepeatMasker"
vx "repeatmasker" $(RepeatMasker  | head -n1 |  sed 's/[A-Za-z\+]//g' | sed 's/  -//') 0 

cat $CFG_PATH | grep "blastdb:" &> /dev/null 

if [ $? -eq "0" ]
then
    ex "blastn"
    vx "blast" $( blastn -version  |head -n 1| awk '{print $2;}' | sed 's/[A-Za-z\+]//g') 1
fi

#cat $CFG_PATH | grep 



if [[ "$@" == *"cluster-config"* ]]; then
    mkdir -p ~/.slurm-logs
    snakemake -s "${SCRIPT_PATH}/Snakefile" -d ${SCRIPT_PATH}  "$@" --cluster "sbatch --parsable -J {cluster.name} -p {cluster.part} -c {cluster.c} --mem {cluster.mem}  -t {cluster.time} --output '$HOME/.slurm-logs/%j.out' --error '$HOME/.slurm-logs/%j.err'";
else
    snakemake -s "${SCRIPT_PATH}/Snakefile" -d ${SCRIPT_PATH}  "$@";
fi

