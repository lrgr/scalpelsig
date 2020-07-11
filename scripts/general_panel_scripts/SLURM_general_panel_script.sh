#!/bin/bash

# TO RUN THIS SCRIPT, do:
# sbatch SLURM_general_panel_script.sh <signature number> <file_tag>

# Lines that begin with #SBATCH specify commands to be used by SLURM for scheduling

#SBATCH --job-name=compute_general_panel                                   # sets the job name
#SBATCH --output out/SCRIPT_OUTS/GENERAL_PANEL_DIRS/GP_SLURM_LOGS/general_panel_log_jobid%j.txt                             # indicates a file to redirect STDOUT to; %j is the jobid
#SBATCH --error out/SCRIPT_OUTS/GENERAL_PANEL_DIRS/GP_SLURM_LOGS/general_panel_error_jobid%j.txt                              # indicates a file to redirect STDERR to; %j is the jobid
#SBATCH --qos=throughput
#SBATCH --time=18:00:00                                         # how long you think your job will take to complete; format=hh:mm:ss
#SBATCH --nodes=1                                               # number of nodes to allocate for your job
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1                                     # request 1 cpu core be reserved per task
#SBATCH --mem 36gb                                               # memory required by job; if unit is not specified MB will be assumed

module add R                                        # run any commands necessary to setup your environment
module use /cbcb/sw/RedHat-7-x86_64/users/jfan03/modules
module add conda
. /cbcb/sw/RedHat-7-x86_64/users/jfan03/local/conda/4.3.24/etc/profile.d/conda.sh
module add R/common/3.5.1
conda activate signature-panel-env

for i in {1..4}; do
    for obj in {1..3}; do
        srun -N 1 -n 1 --mem=35gb --exclusive Rscript find_windows_general_panel.R -s ${1} -t ${2}_it${i} -o ${obj} -w 250 &
    done
done

wait
