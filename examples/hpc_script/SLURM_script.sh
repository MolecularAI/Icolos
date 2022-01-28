#!/bin/bash -l
#SBATCH -N 1
#SBATCH -t 0-02:59:00
#SBATCH -p core
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=2G

source /projects/cc/mai/miniconda3/bin/activate /projects/cc/mai/miniconda3/envs/Icolos
python /projects/cc/mai/Icolos/executor.py -conf /projects/cc/mai/examples/Icolos/MPI_test/workflow_ReSCoSS.json \
       --global_variables "output_dir:<path>/icolos/tests/junk" -debug

