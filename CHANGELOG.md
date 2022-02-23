# Update log

### Version 1.7.1 | 2022-02-23
#### Features
- Added version logging.
#### Internal
- Update of environment specification.
- Improvements to `ESPsim` integration.

### Version 1.7.0 | 2022-02-16
#### Features
- Added Gold docking backend.
- Added OpenFF parametrisation for GROMACS simulations.
#### Internal
- Refactor MD system parametrisation and GROMACS topology handling.
- Refactor step dispatching to allow for cleaner parallelization.
- Various bug fixes and stability improvements.

### Version 1.6.0 | 2022-02-04
#### Internal
- Unit tests fixed.
- Improved stability for `AutoDock Vina` and `RDkit` embedding steps.
- Change nomenclature for parallel step dispatch and Slurm interface
- Improve error handling for PMX steps

### Version 1.5.0 | 2022-02-01
#### Internal
- Refactored gromacs + free energy tools to use ambertools installed in the conda env.
- Restructured repository to make the package `pip` installable.
- Improved unit test coverage for pmx steps.

### Version 1.4.0 | 2022-01-19
#### Features
- Added support for non-equilibrium relative binding free energy calculation with `PMX`.
- Added Glide support for feeding in "in" files from Maestro directly.
- Added AutoDock Vina as docking backend.

#### Internal
- Limited refactoring of support functions.

### Version 1.3.0 | 2021-11-18
#### Features
- Added Iterator mechanism for parallel step execution.
- Pose rescoring my RMSD workflow.
- MMGBSA workflow with GROMACS.

#### Internal
- Improved error logging from subprocesses.
- Improvements to MDpocket workflows.
- Refactored example workflows + added new examples.

### Version 1.2.0 | 2021-09-15
#### Features
- Added MDpocket workflow for pocket identification.
- Expanded scope of GROMACS workflow for improved ligand/cofactor parametrisation.
- Improved FEP+ workflow map construction and analysis.
- Performance optimisation for Turbomole and Prime.
- Added PDBFixer step.
- Added ensemble docking.

#### Internal
- Improved temporary file handling.

### Version 1.1.0 | 2021-06-30
#### Features
- Added `Ligprep` workflow step.
- Added `Glide` workflow step.
- Added run-time global variables.
- Added JSON input type (`REINVENT`-compatible).
- Additional `GROMACS` binaries, and automated ligand parametrisation.
- Added support for Schrodinger's `FEP+` workflow.
- Added support for `OptunaAZ` model building.

#### Bug fixes
- Fixed problems in tabular write-out (no compound names and sometimes lost column order).
- Fixed bug in aggregation (`highest_is_best` parameter was not working properly).
- Fixed instability with step write-out (occurred when no conformers were associated with a compound).
- Fixed bug in the parallelization of `Ligprep`.

#### Internal
- Refactored structure for `Schrodinger` binaries.
- Reworked the write-out functionality.
- Reworked internal file handling.
- Reworked generic data handling.

### Version 1.0.0 | 2021-05-21
#### Features
- Basic functionality (data handling, backend wrapping).
- Various steps implemented (`Turbomole`, `Cosmo`, `OMEGA`, `GROMACS`, ...).

### Bug fixes
- Fixed issues with `Turbomole` execution.
- Enforced GROMACS execution in `tmp_dir`.

### Internal
- Adapted `pydantic` interface.

