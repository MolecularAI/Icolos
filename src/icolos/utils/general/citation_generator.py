import os
from typing import List


from icolos.core.workflow_steps.step import StepBase

def print_citations(steps: List[StepBase], logger = None):
    header = """
=====================================================================
ooooo   .oooooo.     .oooooo.   ooooo          .oooooo.    .oooooo..o 
`888'  d8P'  `Y8b   d8P'  `Y8b  `888'         d8P'  `Y8b  d8P'    `Y8 
 888  888          888      888  888         888      888 Y88bo.      
 888  888          888      888  888         888      888  `"Y8888o.  
 888  888          888      888  888         888      888      `"Y88b 
 888  `88b    ooo  `88b    d88'  888       o `88b    d88' oo     .d8P 
o888o  `Y8bood8P'   `Y8bood8P'  o888ooooood8  `Y8bood8P'  8""88888P'  
=====================================================================\n"""
    writeout_lines = header.split('\n')
    writeout_lines.append("If you publish work using Icolos, please consider citing the following papers, generated based on the workflow steps...\n\n")
    citations = []
    for step in steps:
        if "gmx" in step.type:
            citations.append("\nGROMACS: https://doi.org/10.1016/j.softx.2015.06.001\nH.J.C. Berendsen, D. van der Spoel, and R. van Drunen, “GROMACS: A message-passing parallel molecular dynamics implementation,” Comp. Phys. Comm., 91 43–56 (1995)\nD. van der Spoel, E. Lindahl, B. Hess, G. Groenhof, A.E. Mark, and H.J.C. Berendsen, “GROMACS: Fast, Flexible and Free,” J. Comp. Chem., 26 1701–1718 (2005)\nH. Bekker, H.J.C. Berendsen, E.J. Dijkstra, S. Achterop, R. van Drunen, D. van der Spoel, A. Sijbers, and H. Keegstra et al., “Gromacs: A parallel computer for molecular dynamics simulations”; pp. 252–256 in Physics computing 92. Edited by R.A. de Groot and J. Nadrchal. World Scientific, Singapore, 1993")
        if "pmx" in step.type:
            citations.append("\nPMX: Vytautas Gapsys, Servaas Michielssens, Daniel Seeliger and Bert L. de Groot. Accurate and Rigorous Large Scale Mutation Free Energy Prediction. Angewandte Chemie Int. Ed. 55: 7364-7368(2016)\nVytautas Gapsys, Servaas Michielssens, Daniel Seeliger, and Bert L. de Groot. pmx: Automated protein structure and topology generation for alchemical perturbations. J. Comput. Chem. 36:348-354 (2015\n")
        if "MMPBSA" in step.type:
            citations.append("\ngmx_MMPBSA: Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS. Journal of Chemical Theory and Computation, 2021 17 (10), 6281-6291\n")
        if "turbomole" in step.type:
            citations.append("\nTURBOMOLE: Modular program suite for ab initio quantum-chemical and condensed-matter simulations; J. Chem. Phys. 152, 184107 (2020); https://doi.org/10.1063/5.0004635\n")
        if "espsim" in step.type:
            citations.append("\nESPSim: https://doi.org/10.26434/chemrxiv-2021-sqvv9-v3")
        if "vina" in step.type:
            citations.append("\nAutoDock Vina: Eberhardt, J., Santos-Martins, D., Tillack, A.F., Forli, S. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. Journal of Chemical Information and Modeling.\nTrott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. Journal of computational chemistry, 31(2), 455-461.")
            
    citations = list(set(citations))
    for c in citations:
        writeout_lines.append(c)
        
        
            
    for line in writeout_lines:
        print(line.center(os.get_terminal_size().columns))
    if logger is not None:
        for line in header.split("\n"):
            logger.log(line)
    



    