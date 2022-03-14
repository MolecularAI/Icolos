[![License](https://img.shields.io/badge/License-Apache_2.0-green.svg)](https://opensource.org/licenses/Apache-2.0)
[![PyPI version](https://badge.fury.io/py/icolos.svg)](https://badge.fury.io/py/icolos)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black) 
[![GitHub contributors](https://badgen.net/github/contributors/MolecularAI/Icolos)](https://GitHub.com/MolecularAI/Icolos/graphs/contributors/)
[![Latest tag](https://badgen.net/github/tag/MolecularAI/Icolos)](https://github.com/MolecularAI/Icolos/tag)
[![GitHub forks](https://img.shields.io/github/forks/MolecularAI/Icolos.svg?style=social&label=Fork&maxAge=2592000)](https://GitHub.com/MolecularAI/Icolos/network/)
[![GitHub stars](https://img.shields.io/github/stars/MolecularAI/Icolos.svg?style=social&label=Star&maxAge=2592000)](https://GitHub.com/MolerularAI/Icolos/stargazers/)


# `Icolos`: Workflow manager

The `Icolos` tool is a workflow manager for structure-based workflows in computational chemistry, that abstracts execution logic from implementation as much as possible. Icolos was designed to interface with [REINVENT](https://github.com/MolecularAI/Reinvent), and workflows can be called as a component of the scoring function, or to postprocess results with more expensive methods. Workflows are specified in `JSON` format (see folder `examples`). Currently wrapped are a diverse set of tools and internal steps, including docking, QM and MD capabilities. The pre-print is available [here](https://doi.org/10.26434/chemrxiv-2022-sjcp3).


## Introduction
`Icolos` provides a unified interface to a host of software for common computational chemistry calculations, with built in parallelization,
and straight-forward extensibiltiy to add additional functionality. It was principally developed to handle structural calculations for `REINVENT` jobs, however, workflows can also be run independently.

Workflows are constructed from elementary 'steps', individual blocks which are combined to specify a sequential list of operations, with control of the command-line options provided through step settings, and options to control other aspects of the step's behaviour included in the `additional` block.

For many use cases, one of the template workflows might suit your needs, or need a few tweaks to do what you want. Demonstration notebooks for common workflows are available [here](https://github.com/MolecularAI/IcolosCommunity).

## Initial configuration
You are welcome to clone the repository and use a local version, and in particular if you would like to experiment with the code base and/or contribute features, please get 
in contact with us.

## Installation
After cloning, first install and activate the `icolosprod` `conda` environment:
```
conda env create -f environment_min.yml
conda activate icolosprod
```
Then install the package:
```
pip install -e .
```
This will give you access to the `icolos` entrypoint.

### `ESPsim` installation
The following will install the `ESPsim` package into the environment - this is only required if ligand-based matching using this package is desired.

```
cd ..
git clone https://github.com/hesther/espsim.git
cd espsim
conda activate icolosprod
pip install -e .
```
## Unit testing
Icolos is extensively unit tested, and relies on an external data repo located [here](https://github.com/MolecularAI/IcolosData).  The full test suite takes ~60 mins on a workstation, therefore it is recommended that you execute a subset of unit tests relevant to the workflow you are running.  To execute the full test suite, run something like:
```
pytest -n 8 tests/
```

## Execution
Once a `JSON` is specified, the workflow can be executed like so:

```
conda activate icolosprod
icolos -conf workflow.json
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
icolos -conf workflow.json
```
For GROMACS workflows requiring the GPU partition, you will need to adapt the header accordingly, e.g. like so:

```
#!/bin/bash
#SBATCH -J gmx_cco1_fold_microsecond
#SBATCH -o MygpuJob_out_%j.txt
#SBATCH -e MygpuJob_err_%j.txt
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=4g
#SBATCH -p gpu
#SBATCH --time=12:00:00

```

## Developers
- Christian Margreitter [@cmargreitter](https://github.com/CMargreitter)
- J. Harry Moore [@jharrymoore](https://github.com/jharrymoore)
- Matthias R. Bauer <mattias.r.b@gmail.com>
