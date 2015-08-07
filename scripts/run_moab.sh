#!/bin/bash

# submit via (cmdline overrides script):
# msub -q genacc_q -N latgas_cmdln moab.sh

#queues = genacc_q, backfill, backfill2

#MOAB -N "latgas2"
#MOAB -l nodes=16
#MOAB -j oe
#MOAB -q genacc_q
#MOAB -m abe
#MOAB -l walltime=24:00:00



echo "------------------------------------------------------"
echo "-n 'Job is running on node '; cat $PBS_NODEFILE"
echo "------------------------------------------------------"
echo "MOAB: the queue manager (qsub) is running on $PBS_O_HOST"
echo "MOAB: this job was submitted to the queue: $PBS_O_QUEUE"
echo "MOAB: this job is actually running on the queue: $PBS_QUEUE"
echo "MOAB: working directory is: $PBS_O_WORKDIR"
echo "MOAB: execution mode is: $PBS_ENVIRONMENT"
echo "MOAB: job identifier is: $PBS_JOBID"
echo "MOAB: job name is: $PBS_JOBNAME"
echo "MOAB: node file is: $PBS_NODEFILE"
echo "MOAB: current home directory is: $PBS_O_HOME"
echo "MOAB: PATH = $PBS_O_PATH"
echo "------------------------------------------------------"
echo "-n 'Job is running on node '; cat $PBS_NODEFILE"
echo "------------------------------------------------------"
echo ' '
echo ' '

cd $PBS_O_WORKDIR
pwd
cp $PBS_NODEFILE ./pbs_nodes
echo "NODES: BEGIN"
cat ./pbs_nodes
echo "NODES: END"
date
module load gnu-openmpi
ls -l
mpirun --hostfile ./pbs_nodes ./Latgas2.exe --numwin=4 --numwalk=8 --Elo=-9 --Ehi=-4
date
