
import shutil
from collections import defaultdict
import statistics
from pathlib import Path
import random

import pandas as pd

from evolvex.model_dataclasses import MC_Model
from evolvex.foldx_commands import create_individual_list_foldx_mutations_file, run_foldx_BuildModel, get_complex_stability_ddG, get_chain_group_stability_dG
from evolvex.utils_bio import get_chain_to_sequence_map

large_hydrophobic_residues = 'FILWY'

def create_model(input_PDB_file_path, copy_PDB_file_to_output_dir, mutations_list, output_dir, GLOBALS, output_file_tag=None, PDB_file_tag=None):
    """
    Creates a PDB model using BuildModel, taking as input a PDB file and a list of mutations.
    """
    output_dir.mkdir(exist_ok = True)

    if copy_PDB_file_to_output_dir:
        input_PDB_file_path = shutil.copy(
            src = input_PDB_file_path, 
            dst = output_dir
        )
        input_PDB_file_path = Path(input_PDB_file_path)

    individual_list_foldx_mutations_file_path = create_individual_list_foldx_mutations_file(
        mutant = ','.join(mutations_list), 
        output_dir = output_dir
    )
    
    run_foldx_BuildModel(
        foldx_dir=GLOBALS.foldx_dir, 
        PDB_file_dir=input_PDB_file_path.parent, PDB_file_name=input_PDB_file_path.stem,
        individual_list_foldx_mutations_file_path=individual_list_foldx_mutations_file_path,
        move_neighbors_flag=True,
        vdwDesign=GLOBALS.vdwDesign,
        print_stdout=GLOBALS.print_stdout,
        output_dir=output_dir, output_file_tag=output_file_tag, PDB_file_tag=PDB_file_tag,
    )

    return


def get_acceptable_positions_mut_names_map(all_mutations_summary_df, PDB_name):
    acceptable_mutations_map = defaultdict(list)
    for position, position_df in all_mutations_summary_df.groupby('position'):
        binding_ddG_variance = statistics.variance(position_df.binding_ddG.values)
        antibody_stability_ddG_variance = statistics.variance(position_df.antibody_stability_ddG.values)
        mean_binding_and_stability_variance = (binding_ddG_variance + antibody_stability_ddG_variance) / 2
        for mut_name, row in position_df.iterrows():
            mutant_residue = mut_name[-1]
            # This filters out mutations that are too destabilizing either in terms of binding or stability
            if row.complex_stability_ddG > 1 or row.binding_ddG > 1 or row.antibody_stability_ddG > 2:
                continue

            # We don't want mutations to large hydrophobics when mutations at that position don't seem to change binding or stability a lot (ddG values ~ 0)
            if mean_binding_and_stability_variance < 0.1 and mutant_residue in large_hydrophobic_residues:
                continue
            
            acceptable_mutations_map[position].append(mut_name)

        # If no mutation passes the filters, allow all mutations at that position and let the search algorithm find what's best
        if not position in acceptable_mutations_map:
            print(f'No acceptable mutations found for {position = } in {PDB_name = }. Allowing all mutations at that position during the search. Probably needs manual inspection.', flush=True)
            acceptable_mutations_map[position].extend(position_df.index.values)

    return acceptable_mutations_map

def get_hotspot_positions_mut_names_map(all_mutations_summary_df):
    hotspot_mutations_map = defaultdict(list)
    for mut_name, row in all_mutations_summary_df.iterrows():
        original_residue = row.original_residue
        mutant_residue = mut_name[-1]
        position = mut_name[2:-1]

        # Hotspot mutations need to strongly improve binding affinity and not be too destabilizing for the antibody
        if row.binding_ddG < -1.5 and row.antibody_stability_ddG < 1:
            hotspot_mutations_map[position].append(mut_name)

        # If mutating from Alanine to the original residue is highly stabilizing for the antibody, then we consider it as a hotspot mutation, even if it maybe
        # doesn't contribute to binding
        elif mutant_residue == original_residue and row.antibody_stability_ddG < -2:
            hotspot_mutations_map[position].append(mut_name)

    return hotspot_mutations_map

def get_allowed_mutations_per_position_maps(PDB_name, all_mutations_summary_file_path):
    """
    Returns two dictionaries:
    1) A dict where keys are positions and values are mutation names (i.e {'53':['A53R', 'A53K', ...]})
    2) A dict where values are amino acid mutations (i.e {'53':['R', 'K', ...]}) 
    """
    all_mutations_summary_df = pd.read_csv(all_mutations_summary_file_path, header=0, index_col=0, dtype = {'position':str})

    acceptable_positions_mut_names_map = get_acceptable_positions_mut_names_map(all_mutations_summary_df, PDB_name)
    hotspot_positions_mut_names_map = get_hotspot_positions_mut_names_map(all_mutations_summary_df)

    allowed_mut_names_per_position_map = {}
    allowed_AA_per_position_map = {}
    for position, position_df in all_mutations_summary_df.groupby('position'):
        # In positions with an antibody folding hotspot mutation, only consider hotspot mutations, otherwise allow both hotspot and acceptable mutations.
        if any(position_df.antibody_stability_ddG < -2):
            allowed_mut_names = hotspot_positions_mut_names_map[position]
            allowed_AA_mutations = {mut_name[-1] for mut_name in allowed_mut_names}
        else:
            allowed_mut_names = hotspot_positions_mut_names_map[position] + acceptable_positions_mut_names_map[position]
            allowed_AA_mutations = {mut_name[-1] for mut_name in allowed_mut_names}

        if allowed_mut_names:
            allowed_mut_names_per_position_map[position] = allowed_mut_names
            allowed_AA_per_position_map[position] = allowed_AA_mutations
        
    return allowed_mut_names_per_position_map, allowed_AA_per_position_map 


def clean_up_model_dir(model_dir, PDB_file_name_to_keep_as_model):
    for file in model_dir.iterdir():
        if file.name != PDB_file_name_to_keep_as_model:
            file.unlink()

    PDB_file_path_to_keep_as_model = model_dir / PDB_file_name_to_keep_as_model
    if PDB_file_name_to_keep_as_model != 'model.pdb':
        PDB_file_path_to_keep_as_model.rename(PDB_file_path_to_keep_as_model.with_name('model.pdb'))
    return

def get_random_mutations_list_for_initial_population(allowed_mut_names_per_position_map):
    random_mutations_list = []
    for position, mut_names_list in allowed_mut_names_per_position_map.items():
        # Skip positions where the wildtype residue is not Alanine, as this means the position was marked as MakeAla = "N"
        if all(mut_name[0] != 'A' for mut_name in mut_names_list):
            continue

        random_mut_name = random.choice(mut_names_list)
        random_mutations_list.append(random_mut_name)
    
    random.shuffle(random_mutations_list)
    return random_mutations_list

def generate_random_model(
        PDB_name, foldx_Alanine_mutant_PDB_file_path, allowed_mut_names_per_position_map, allowed_AA_mutations_per_position_map,  antibody_stability_dG_original_wildtype,
        antibody_seq_map_original_wildtype, model_dir, GLOBALS
    ):
    model, n_tries = None, 0
    while model == None:
        random_mutations_list = get_random_mutations_list_for_initial_population(allowed_mut_names_per_position_map)

        create_model(
            input_PDB_file_path = foldx_Alanine_mutant_PDB_file_path, 
            copy_PDB_file_to_output_dir = True, 
            mutations_list = random_mutations_list, 
            output_dir = model_dir, 
            GLOBALS = GLOBALS
        )

        complex_stability_ddG = get_complex_stability_ddG(model_dir / f'Dif_{PDB_name}_1_Alanine_mutant.fxout')
        if complex_stability_ddG < 0.5 or n_tries == 5:
            full_residue_IDs_list = [f'{mut_name[-1]}{mut_name[1:-1]}' for mut_name in random_mutations_list]
            model = MC_Model(
                model_dir = model_dir,
                full_residue_IDs_list = full_residue_IDs_list,
                backbone_PDB_file_name = PDB_name,
                antibody_stability_dG_original_wildtype = antibody_stability_dG_original_wildtype,
                antibody_seq_map_original_wildtype = antibody_seq_map_original_wildtype,
                allowed_AA_mutations_per_position_map = allowed_AA_mutations_per_position_map,
            )
        else:
            shutil.rmtree(model_dir) # Empty the directory and try again with a new random list of mutations

        n_tries += 1

    clean_up_model_dir(model_dir = model_dir, PDB_file_name_to_keep_as_model = f'{PDB_name}_1_Alanine_mutant_1.pdb')
    return model


def generate_initial_models(parallel_executor, evolvex_working_dir, backbone_PDB_files_paths, GLOBALS):
    futures = []
    for PDB_file_path in backbone_PDB_files_paths:
        PDB_dir = evolvex_working_dir / PDB_file_path.stem
        PDB_name = PDB_dir.name

        foldx_Alanine_mutant_PDB_file_path = PDB_dir / f'{PDB_name}_1_Alanine_mutant.pdb'

        search_output_dir = PDB_dir / 'search_results'; search_output_dir.mkdir(exist_ok = True)

        all_mutations_summary_file_path = PDB_dir / 'hotspot_mutants' / 'all_mutations_summary.csv'
        if not all_mutations_summary_file_path.exists():
            print(f'Could not find the all_mutations_summary.csv file for {PDB_name = }, this should not happen ! Skipping PDB backbone for search.', flush=True)
            continue

        allowed_mut_names_per_position_map, allowed_AA_mutations_per_position_map = get_allowed_mutations_per_position_maps(
            PDB_name = PDB_name, all_mutations_summary_file_path = all_mutations_summary_file_path
        )

        antibody_stability_dG_original_wildtype = get_chain_group_stability_dG(indiv_file_path = PDB_dir / 'Indiv_energies_original_wildtype_AC.fxout', chain_group_name = GLOBALS.antibody_chains)
        antibody_seq_map_original_wildtype = get_chain_to_sequence_map(PDB_file_path = foldx_Alanine_mutant_PDB_file_path, chain_subset = GLOBALS.antibody_chains)

        for ith_model in range(GLOBALS.population_size):
            model_dir = search_output_dir / str(ith_model); model_dir.mkdir(exist_ok=True)
            future = parallel_executor.submit(
                generate_random_model, PDB_name, foldx_Alanine_mutant_PDB_file_path, allowed_mut_names_per_position_map, allowed_AA_mutations_per_position_map,  
                antibody_stability_dG_original_wildtype, antibody_seq_map_original_wildtype, model_dir, GLOBALS
            )
            futures.append(future)

    initial_models_population = parallel_executor.gather(futures)
    
    return initial_models_population