[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black) 


# `Icolos`: Workflow manager

The `Icolos` tool is a workflow manager, that aims at separating execution logic from actual implementation as much as
possible. Workflows are specified in `JSON` files (see folder `examples`), linking steps together. Currently wrapped are
a diverse set of tools and internal steps, including QM and MD software.


## Introduction
`Icolos` provides a single, unified interface to a host of software for common computational chemistry calculations, with built in parallelization,
and straight-forward extensibiltiy to add additional functionality. It was principally developed to handle structural calculations for `REINVENT` jobs, however, various workflows have also been used as stand-alone pipelines.

Workflows are constructed from elementary 'steps', individual blocks which specify a sequential list of operations (normally corresponding to a single program being executed),
with control of the command-line options provided through step settings, and options to control other aspects of the step's behaviour included in the `additional` block.

For many use cases, one of the template workflows might suit your needs, or need a few tweaks to do what you want. The JSONs in the example folder are less regularly updated
with new features and are mostly used for integration testing.

## Initial configuration
You are welcome to clone the repository and use a local version, and in particular if you would like to experiment with the code base and/or contribute features, please get 
in contact with us.

## Installation
After cloning, first install the `icolosprod` `conda` environment:
```
conda create -f environment_min.yml
```

## Execution
Once a `JSON` is specified, the workflow can be executed like so:

```
conda activate icolosprod
python executor.py -conf workflow.json
```

## `SLURM` Execution
Once specified, a workflow can be called like this in a `bash` script:

```
#!/bin/bash -l
#SBATCH -N 1
#SBATCH -t 0-02:59:00
#SBATCH -p core
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=2G

source /<conda_path>/miniconda3/bin/activate /<conda_path>/minconda3/envs/icolosprod
python /<icolos_path>/Icolos/executor.py -conf workflow.json
```
For GROMACS workflows requiring the GPU partition, you will need to adapt the header accordingly, e.g. like so:

```
#!/bin/bash
#SBATCH -J gmx_cco1_fold_microsecond
#SBATCH -o MygpuJob_out_%j.txt
#SBATCH -e MygpuJob_err_%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:4
#SBATCH --gres-flags=enforce-binding
#SBATCH --mem-per-cpu=4g
#SBATCH -p gpu
#SBATCH --time=370:00:00

```

## Developers
- Christian Margreitter <christian.margreitter@astrazeneca.com>
- J. Harry Moore <harry.moore@astrazeneca.com>
- Matthias R. Bauer <mattias.r.b@gmail.com>
