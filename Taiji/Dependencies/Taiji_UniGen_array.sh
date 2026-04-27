#!/bin/bash
#SBATCH -J Taiji_UniGen_Array                    # '##' means comment, 1 '#SBATCH' can be treated as the setting parameters
###SBATCH -N 1
###SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH -t 100:00:00                # set time limit
#SBATCH -a 1-2%1  # you have 138 jobs to run, 15 jobs can be submitted at the same time.
##SBATCH --nodelist=compute-1
input_fp=$1

PARAM1=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" ${input_fp} | awk '{print $1}')
echo $PARAM1
input_dir=$(dirname "$PARAM1")
echo ${input_dir}

# change_one=${input_dir}/Input
# echo $change_one
change_two=${input_dir}/../../Output/Partial/$(basename $PARAM1 _config.yml)_output
echo $change_two
#change_two=${change_one/input/output}
cd $change_two

ml load taiji/1.3.0
wait 

taiji run --config $PARAM1 #--cloud
