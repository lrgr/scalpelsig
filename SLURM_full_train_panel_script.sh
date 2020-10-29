#!/bin/bash

# TO RUN THIS SCRIPT, do:
# sbatch SLURM_full_train_panel_script.sh <file_tag> 

# Lines that begin with #SBATCH specify commands to be used by SLURM for scheduling

#SBATCH --job-name=compute_ft_panel                                   # sets the job name
#SBATCH --output out/SCRIPT_OUTS/SLURM_LOGS/ft_panel_log_jobid%j.txt                             # indicates a file to redirect STDOUT to; %j is the jobid
#SBATCH --error out/SCRIPT_OUTS/SLURM_LOGS/ft_panel_error_jobid%j.txt                              # indicates a file to redirect STDERR to; %j is the jobid
#SBATCH --qos=throughput
#SBATCH --time=18:00:00                                         # how long you think your job will take to complete; format=hh:mm:ss
#SBATCH --nodes=1                                               # number of nodes to allocate for your job
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=1                                     # request 1 cpu core be reserved per task
#SBATCH --mem 36gb                                               # memory required by job; if unit is not specified MB will be assumed

module add R                                        # run any commands necessary to setup your environment
module use /cbcb/sw/RedHat-7-x86_64/users/jfan03/modules
module add conda
. /cbcb/sw/RedHat-7-x86_64/users/jfan03/local/conda/4.3.24/etc/profile.d/conda.sh
module add R/common/3.5.1
conda activate signature-panel-env

for i in 2 3 8 13 18; do
        srun -N 1 -n 1 --mem=35gb --exclusive Rscript full_panel_script.R -s ${i} -t ${1} -o 2 -w 250 &
done

wait
