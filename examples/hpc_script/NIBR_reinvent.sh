#!/bin/bash -l
#SBATCH -N 1
#SBATCH -t 0-02:59:00
#SBATCH -p core
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=2G

source /projects/cc/mai/miniconda3/bin/activate /projects/cc/mai/miniconda3/envs/icolosprod
icolos -conf /projects/cc/mai/material/Icolos/templates/NIBR/12-06-21nibr.json -debug \
       --global_variables "entrypoint_dir:<path>/icolos, input_path_json:{entrypoint_dir}/tests/data/reinvent/small_input.json, output_path_json:{entrypoint_dir}/tests/junk/nibr_reinvent.json"

