import os
from typing import List
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.step_enums import StepBaseEnum

_SBE = StepBaseEnum


class ConsoleColours:
    HEADER = "\033[95m"
    BLUE = "\033[94m"
    CYAN = "\033[96m"
    GREEN = "\033[92m"
    RED = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    ORANGE = "\033[0;33m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    BLINKING = "\33[5m"


def add_citation(
    step_type: str,
    workflow_step_types: List[str],
    citations: List[str],
    citation_string: str,
):
    if any(
        [
            True
            for w_step_type in workflow_step_types
            if step_type.upper() in w_step_type
        ]
    ) or step_type == "default":
        citations.append(citation_string)


def print_citations(steps: List[StepBase], logger=None):
    header = f"""{ConsoleColours.GREEN}
=====================================================================
ooooo   .oooooo.     .oooooo.   ooooo          .oooooo.    .oooooo..o 
`888'  d8P'  `Y8b   d8P'  `Y8b  `888'         d8P'  `Y8b  d8P'    `Y8 
 888  888          888      888  888         888      888 Y88bo.      
 888  888          888      888  888         888      888  `"Y8888o.  
 888  888          888      888  888         888      888      `"Y88b 
 888  `88b    ooo  `88b    d88'  888       o `88b    d88' oo     .d8P 
o888o  `Y8bood8P'   `Y8bood8P'  o888ooooood8  `Y8bood8P'  8""88888P'  
====================================================================={ConsoleColours.ENDC}\n"""
    writeout_lines = header.split("\n")
    writeout_lines.append(
        "If you publish work using Icolos, please consider citing the following papers, based on the workflow's steps...\n\n"
    )
    citations = []
    step_types = [step.type.upper() for step in steps]

    # add general Icolos citation
    add_citation(
        "default",
        step_types,
        citations,
        "\nIcolos: J. Harry Moore, Matthias R. Bauer, Jeff Guo, Atanas Patronov, Ola Engkvist and Christian Margreitter. Icolos: A workflow manager for structure based post-processing of de novo generated small molecules. https://doi.org/10.26434/chemrxiv-2022-sjcp3, https://github.com/MolecularAI/Icolos (2022)\n",
    )

    # add individual citations based on the step types used in this workflow (can be empty, of course)
    add_citation(
        "mdrun",
        step_types,
        citations,
        "\nGROMACS: https://doi.org/10.1016/j.softx.2015.06.001\nH.J.C. Berendsen, D. van der Spoel, and R. van Drunen, “GROMACS: A message-passing parallel molecular dynamics implementation,” Comp. Phys. Comm., 91 43–56 (1995)\nD. van der Spoel, E. Lindahl, B. Hess, G. Groenhof, A.E. Mark, and H.J.C. Berendsen, “GROMACS: Fast, Flexible and Free,” J. Comp. Chem., 26 1701–1718 (2005)\nH. Bekker, H.J.C. Berendsen, E.J. Dijkstra, S. Achterop, R. van Drunen, D. van der Spoel, A. Sijbers, and H. Keegstra et al., “Gromacs: A parallel computer for molecular dynamics simulations”; pp. 252–256 in Physics computing 92. Edited by R.A. de Groot and J. Nadrchal. World Scientific, Singapore, 1993",
    )
    add_citation(
        "pmx_",
        step_types,
        citations,
        "\nPMX: Vytautas Gapsys, Servaas Michielssens, Daniel Seeliger and Bert L. de Groot. Accurate and Rigorous Large Scale Mutation Free Energy Prediction. Angewandte Chemie Int. Ed. 55: 7364-7368(2016)\nVytautas Gapsys, Servaas Michielssens, Daniel Seeliger, and Bert L. de Groot. pmx: Automated protein structure and topology generation for alchemical perturbations. J. Comput. Chem. 36:348-354 (2015\n",
    )
    add_citation(
        "gmx_MMPBSA",
        step_types,
        citations,
        "\ngmx_MMPBSA: Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS. Journal of Chemical Theory and Computation, 2021 17 (10), 6281-6291\n",
    )
    add_citation(
        "fep",
        step_types,
        citations,
        "\nFEP+: R. Abel, L. Wang, E.D. Harder, B.J. Berne, R.A. Friesner. Advancing Drug Discovery through Enhanced Free Energy Calculations. Acc. Chem. Res., 2017 50 (7), 1625-1632\n",
    )
    add_citation(
        _SBE.STEP_KALLISTO,
        step_types,
        citations,
        "\nKallisto: Caldeweyher, E., Kallisto: A command-line interface to simplify computational modelling and the generation of atomic features. Journal of Open Source Software, 6(60), 3050, https://doi.org/10.21105/joss.03050. 2021\n",
    )
    add_citation(
        _SBE.STEP_JAZZY,
        step_types,
        citations,
        "\nJazzy: Caldeweyher, E., Ghiandoni, G.M., https://github.com/AstraZeneca/jazzy <CITATION MISSING>\n",
    )
    add_citation(
        _SBE.STEP_TURBOMOLE,
        step_types,
        citations,
        "\nTURBOMOLE: Modular program suite for ab initio quantum-chemical and condensed-matter simulations; J. Chem. Phys. 152, 184107 (2020); https://doi.org/10.1063/5.0004635\n",
    )
    add_citation(
        _SBE.STEP_ESP_SIM,
        step_types,
        citations,
        "\nESPSim: https://doi.org/10.26434/chemrxiv-2021-sqvv9-v3\n",
    )
    add_citation(
        _SBE.STEP_AUTODOCKVINA_DOCKING,
        step_types,
        citations,
        "\nAutoDock Vina: Eberhardt, J., Santos-Martins, D., Tillack, A.F., Forli, S. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. Journal of Chemical Information and Modeling.\nTrott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. Journal of computational chemistry, 31(2), 455-461.\n",
    )
    add_citation(
        _SBE.STEP_XTB,
        step_types,
        citations,
        "\nGFN2-xTB: C. Bannwarth, S. Ehlert, S. Grimme. GFN2-xTB—An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions.\nJ. Chem. Theory Comput. 2019, 15 (3), 1652–1671 DOI: 10.1021/acs.jctc.8b01176\n",
    )
    add_citation(
        _SBE.STEP_OMEGA,
        step_types,
        citations,
        "\nOMEGA: P.C.D Hawkins, A.G. Skillman, G.L. Warren, B.A. Ellingson, M.T. Stahl. Conformer Generation with OMEGA: Algorithm and Validation Using High Quality Structures from the Protein Databank and the Cambridge Structural Database. J. Chem. Inf. Model. 2010, 50, 572-584.\n",
    )
    add_citation(
        _SBE.STEP_COSMO,
        step_types,
        citations,
        "\nCOSMOtherm, Version C3.0, Release 17.01; COSMOlogic GmbH & Co. KG, http://www.cosmologic.de\n",
    )
    add_citation(
        _SBE.STEP_SHAEP,
        step_types,
        citations,
        "\nShaEP: M.J. Vainio, J.S. Puranen, M.S. Johnson. ShaEP: Molecular Overlay Based on Shape and Electrostatic Potential. J. Chem. Inf. Model. 49, 492-502 (2009).\n",
    )
    add_citation(
        _SBE.STEP_PANTHER,
        step_types,
        citations,
        "\nS.P. Niinivehmas, K. Salokas, S. Lätti, H. Raunio, O.T. Pentikäinen. Ultrafast protein structure-based virtual screening with Panther. J. Comput. Aided. Mol. Des. 29, 989–1006. doi: 10.1007/s10822-015-9870-3 (2015).\n",
    )
    add_citation(
        _SBE.STEP_MDPOCKET,
        step_types,
        citations,
        "\nMDpocket: P. Schmidtke, A. Bidon-Chanal, F.J. Luque, X. Barril. MDpocket: open-source cavity detection and characterization on molecular dynamics trajectories. Bioinformatics, Volume 27, Issue 23, Pages 3276–3285. https://doi.org/10.1093/bioinformatics/btr550 (2011)\n",
    )
    add_citation(
        _SBE.STEP_GOLD_DOCKING,
        step_types,
        citations,
        "\nG. Jones, P. Willett, R.C. Glen, A.R. Leach, R. Taylor. Development and Validation of a Genetic Algorithm for Flexible Docking. J. Mol. Biol., 267, 727-748. doi: 10.1006/jmbi.1996.0897 (1997).\n",
    )
    add_citation(
        _SBE.STEP_CREST,
        step_types,
        citations,
        "\nCREST: P. Pracht, F. Bohle, S. Grimme. Phys. Chem. Chem. Phys., 22, 7169-7192. doi: 10.1039/C9CP06869D (2020).\n",
    )
    add_citation(
        _SBE.STEP_GLIDE,
        step_types,
        citations,
        "\nGLIDE: R.A. Friesner, R.B. Murphy, M.P. Repasky, L.L. Frye, J.R. Greenwood, T.A. Halgren, P.C. Sanschagrin, D.T. Mainz. Extra Precision Glide: Docking and Scoring Incorporating a Model of Hydrophobic Enclosure for Protein-Ligand Complexes. J. Med. Chem., 49, 6177–6196 (2006).\n",
    )
    add_citation(
        _SBE.STEP_LIGPREP,
        step_types,
        citations,
        "\nSchrödinger Release 2021-4: LigPrep, Schrödinger, LLC, New York, NY (2021).\n",
    )
    add_citation(
        _SBE.STEP_MACROMODEL,
        step_types,
        citations,
        "\nSchrödinger Release 2021-4: MacroModel, Schrödinger, LLC, New York, NY (2021).\n",
    )
    add_citation(
        _SBE.STEP_PRIME,
        step_types,
        citations,
        "\nM.P. Jacobson, D.L. Pincus, C.S. Rapp, T.J.F. Day, B. Honig, D.E. Shaw, R.A. Friesner. A Hierarchical Approach to All-Atom Protein Loop Prediction. Proteins: Structure, Function and Bioinformatics, 55, 351-367 (2004).\n",
    )

    citations = list(set(citations))
    for c in citations:
        writeout_lines.append(c)
    writeout_lines.append("\n")
    writeout_lines.append(
        "=====================================================================\n"
    )

    for line in writeout_lines:
        try:
            width = os.get_terminal_size().columns
        except:
            width = 150

        print(line.center(width))
