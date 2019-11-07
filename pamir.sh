#!/bin/bash
set -e
SCRIPT_PATH="$(dirname `which $0`)/pamir"

#command -v conda && echo "Found conda, skipping checks :)" >&2 || { NO_CONDA=1; }

VERS_CHECKER="./pamir/scripts/version_check.py"

   
#Checking if tools exist
command -v snakemake && echo "... Exists :)" >&2 || { echo snakemake missing! >&2; FAILED=1; }
command -v samtools  && echo "... Exists :)" >&2 || { echo samtools missing!  >&2; FAILED=1; }
command -v blastn    && echo "... Exists :)" >&2 || { echo blastn   missing!  >&2; FAILED=1; }
command -v bedtools  && echo "... Exists :)" >&2 || { echo bedtools missing!  >&2; FAILED=1; }
command -v mrsfast   && echo "... Exists :)" >&2 || { echo mrsfast  missing!  >&2; FAILED=1; } 
#command -v minia   && echo "... Exists :)" >&2 || { echo minia  missing!  >&2; FAILED=1; } 
command -v bwa   && echo "... Exists :)" >&2 || { echo bwa  missing!  >&2; FAILED=1; } 

#Checking Versions are ok

if [ -z $FAILED ]; then 
    BEDTOOLS_VERSION=$( bedtools  --version |head -n 1| awk '{print $2;}' | cut -d"-" -f1 | sed 's/[A-Za-z\+]//g');
    SAMTOOLS_VERSION=$( samtools  --version |head -n 1| awk '{print $2;}' | sed 's/[A-Za-z\+]//g');
    MRSFAST_VERSION=$( mrsfast   --version |head -n 1| awk '{print $2;}' | sed 's/[A-Za-z\+]//g');
    SNAKEMAKE_VERSION=$( snakemake --version |head -n 1| sed 's/[A-Za-z\+]//g');
    BLASTN_VERSION=$( blastn    -version  |head -n 1| awk '{print $2;}' | sed 's/[A-Za-z\+]//g');

#    MINIA_VERSION=$(minia --version | head -n1 | awk '{print $3;}');
    #TODO Check BWA version



    python $VERS_CHECKER gte 5.3.0 $SNAKEMAKE_VERSION && echo "Snakemake version...OK" >&2 || { echo "Snakemake version...should be >= 5.3.0" >&2; FAILED=1; };
    python $VERS_CHECKER gte 1.9 $SAMTOOLS_VERSION && echo "Samtools version...OK" >&2 || { echo "Samtools version...should be >= 1.9" >&2; FAILED=1; };
    python $VERS_CHECKER gte 3.4.0 $MRSFAST_VERSION && echo "MrsFAST version...OK" >&2 || { echo "MrsFAST version...should be >= 3.4.0" >&2; FAILED=1; };
    python $VERS_CHECKER gte 2.26.0 $BEDTOOLS_VERSION && echo "Bedtools version...OK" >&2 || { echo "Bedtools version...should be >= 2.26.0" >&2; FAILED=1; };
    python $VERS_CHECKER gte 2.9.0 $BLASTN_VERSION && echo "Snakemake version...OK" >&2 || { echo "Blastn version...should be >= 2.9.0+" >&2; FAILED=1; };
#    python $VERS_CHECKER gte 3.2.0 $MINIA_VERSION && echo "Minia version...OK" >&2 || { echo "Minia version...should be >= 3.2.0" >&2; FAILED=1; };

fi


if [ -z $FAILED ]; then 
    snakemake -s "${SCRIPT_PATH}/Snakefile" -d ${SCRIPT_PATH}  "$@";
fi
