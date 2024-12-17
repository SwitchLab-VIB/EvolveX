from pathlib import Path
from dataclasses import dataclass

@dataclass
class MC_Model:
    model_dir: Path 
    full_residue_IDs_list: list[str] # Full residue IDs of the mutable positions
    
    # Constants needed for the MC
    backbone_PDB_file_name: str
    antibody_stability_dG_original_wildtype: float
    antibody_seq_map_original_wildtype: dict
    allowed_AA_mutations_per_position_map: dict

@dataclass
class GA_Model:
    model_dir: Path
    full_residue_IDs_list: list[str]

    # Parameters used to generate a Model's PDB file using BuildModel with only 1 mutation based on the parent instead of starting back from the Alanine PDB every time, which would be slow
    parent_model_dir: Path
    mutations_to_generate_PDB: list[str]