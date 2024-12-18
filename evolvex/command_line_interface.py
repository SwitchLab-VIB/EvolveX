import argparse
from pathlib import Path

import argparse
from types import SimpleNamespace
import warnings
import yaml
from pathlib import Path

def read_and_validate_config_file(file_path):
    with open(file_path, 'rt') as f:
        GLOBALS = yaml.safe_load(f)
        GLOBALS = SimpleNamespace(**GLOBALS) # SimpleNamespace enables the use of dot notation to access values instead of dict brackets

    GLOBALS.working_dir = Path(GLOBALS.working_dir)
    GLOBALS.foldx_dir = Path(GLOBALS.foldx_dir)
    GLOBALS.Backbones_dir = Path(GLOBALS.Backbones_dir)
    GLOBALS.PositionsToExplore_file_path = Path(GLOBALS.PositionsToExplore_file_path)
    if not (GLOBALS.working_dir.is_dir()):
        GLOBALS.working_dir.mkdir(parents=True)

    if not (GLOBALS.foldx_dir.is_dir() and (GLOBALS.foldx_dir / 'foldx').exists()):
        raise ValueError("The foldx_dir must contain an executable file named 'foldx'.")

    if GLOBALS.search_algorithm not in ('systematic', 'GA'):
        raise ValueError("The search_algorithm must be one of 'GA' or 'systematic'.")
    
    if GLOBALS.search_algorithm == 'GA':
        if GLOBALS.population_size % 2 != 0:
            raise ValueError("Population_size must be an even number.")
        if GLOBALS.population_size < 50:
            warnings.warn("The population_size for each PDB backbone is < 50, which is low.")

    if GLOBALS.compute_env == 'SLURM':
        SLURM_parameters = ('account_name', 'cluster_name', 'cluster_partition')
        if not all(hasattr(GLOBALS, parameter) for parameter in SLURM_parameters):
            raise ValueError("The SLURM compute_env requires the following parameters: {SLURM_parameters}")

    return GLOBALS

def command_line_interface():
    parser = argparse.ArgumentParser(
        description='EvolveX',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        'config', type=Path,
        help="Path to YAML configuration file. See the evolvex_config_example.yaml file for the list of available parameters, and the README for an exaplanation of each parameter."
    )
    
    args = parser.parse_args()

    GLOBALS = read_and_validate_config_file(file_path = args.config)
    
    return GLOBALS
