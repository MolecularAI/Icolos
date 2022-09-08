PMX 
===

Icolos features a series of steps for automating non-equilibrium binding free energy calculations using Gromacs and PMX

All PMX steps inherit from StepPMXBase:

.. currentmodule:: icolos.core.workflow_steps.pmx

.. autosummary::
   :toctree: generated

	base.StepPMXBase
	setup_workpath.StepPMXSetup
	atomMapping.StepPMXatomMapping
	ligandHybrid.StepPMXligandHybrid
	box_water_ions.StepPMXBoxWaterIons
	assemble_systems.StepPMXAssembleSystems
	prepare_simulations.StepPMXPrepareSimulations
	prepare_transitions.StepPMXPrepareTransitions
	run_simulations.StepPMXRunSimulations
	run_analysis.StepPMXRunAnalysis
	
	

   