from typing import List
from pydantic import BaseModel


class ConfigProteinData(BaseModel):
    protein_datafile: str


class ConfigFitnessFunctionSettings(BaseModel):
    initial_virtual_pt_match_max: int = 3
    relative_ligand_energy: int = 1
    gold_fitfunc_path: str = "goldscore"
    start_vdw_linear_cutoff: float = 6
    score_param_file: str = "DEFAULT"


class ConfigSaveOptions(BaseModel):
    save_score_in_file: int = 1
    save_protein_torsions: int = 1


class ConfigCovalentBonding(BaseModel):
    # TODO: extend to support covalent docking
    covalent: int = 0


class ConfigConstraints(BaseModel):
    # TODO: extend to support constrained docking
    force_constraints: int = 0


class ConfigTermination(BaseModel):
    early_termination: int = 1
    n_top_solutions: int = 3
    rms_tolerance: float = 1.5


class ConfigFlags(BaseModel):
    internal_ligand_h_bonds: int = 0
    flip_free_corners: int = 0
    match_ring_templates: int = 0
    flip_amide_bonds: int = 0
    flip_planar_n: str = "1 flip_ring_NRR flip_ring_NHR"
    flip_pyramidal_n: str = 0
    rotate_carboxylic_oh: str = "flip"
    use_tordist: int = 1
    postprocess_bonds: int = 1
    rotatable_bond_override_file: str = "DEFAULT"
    solvate_all: int = 1


class ConfigDataFiles(BaseModel):
    ligand_data_file: List[str] = []
    param_file: str = "DEFAULT"
    set_ligand_atom_types: int = 1
    set_protein_atom_types: int = 0
    directory: str = "."
    tordist_file: str = "DEFAULT"
    make_subdirs: int = 0
    save_lone_pairs: int = 1
    fit_points_file: str = "fit_pts.mol2"
    read_fitpts: str = 0


class ConfigFloodFill(BaseModel):
    radius: float = 10
    origin: str = "0 0 0"
    do_cavity: int = 1
    floodfill_atom_no: int = 0
    cavity_file: str
    floodfill_center: str = "cavity_from_ligand 7 atoms"


class ConfigGeneticOperators(BaseModel):
    pt_crosswt: str = "auto"
    allele_mutatewt: str = "auto"
    migratewt: str = "auto"


class ConfigPopulation(BaseModel):
    # Note that pydantic accepts also numbers when strings are specified (unless "strict" mode is utilized)
    popsiz: str = "auto"
    select_pressure: str = "auto"
    n_islands: str = "auto"
    maxops: str = "auto"
    niche_siz: str = "auto"


class ConfigAutomaticSettings(BaseModel):
    autoscale: float = (
        1.0  # between 0 (off) and 5.0 (the larger the slower but more sampling)
    )
    autoscale_nops_max: int = (
        100000  # maximum value for autoscale operations (0 means off)
    )
    autoscale_nops_min: int = (
        100  # minimum value for autoscale operations (0 means off)
    )
