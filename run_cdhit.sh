#!/bin/bash
#$ -cwd
#$ -e COMMANDLOGS/
#$ -o COMMANDLOGS/
#$ -q short.qc

# If there's an error, fail the whole script
set -e -o pipefail

# Set permissions so any user in the group can 
# read/write what it's created by the script
umask 002

echo "********************************************************"
# Load software modules
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load networkx/2.2-foss-2019b-Python-2.7.16
echo "Loaded networkx/2.2-foss-2019b-Python-2.7.16 Module"

export PATH=/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/cd-hit-v4.5.7-2011-12-16:$PATH


cdr3file=$1

### Set up parameters 
ID=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$1}" $cdr3file) 
# All Samples for array job
OUTPUT=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$2}" $cdr3file) 

CMD="cd-hit -i ${ID} -d 0 -o ${OUTPUT} -l 5 -n 5 -c 0.85 -S 3 -g 1"
echo "Command : ${CMD}"
echo "********************************************************"
echo
eval "${CMD}"

echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0
