from collections import defaultdict
import itertools
import math
import random

import pandas as pd
from dask.distributed import as_completed

from evolvex.dask_parallel import wait_and_remove
from evolvex.foldx_commands import (get_binding_dG,
                          get_chain_group_intraclash_score, get_complex_stability_dG,
                          get_chain_group_stability_dG,
                          run_foldx_AnalyseComplex,
                          run_foldx_Stability, get_all_other_interaction_file_info)
from evolvex.model_generation import clean_up_model_dir, create_model
from evolvex.model_generation import get_acceptable_positions_mut_names_map, get_hotspot_positions_mut_names_map
from evolvex.utils import save_compressed_PDB_file

from evolvex.utils import NDIGIS_ROUNDING

paratope_AA =         [ 'A', 'C',  'D',  'E', 'F',  'G',  'H',  'I',  'K',  'L', 'M',  'N', 'P', 'Q',  'R',   'S',  'T', 'V', 'W',  'Y']
paratope_AA_weights = [2.95, 0.1, 6.75, 2.85, 3.9, 7.75, 2.95, 2.75, 2.55, 3.55, 0.7, 7.75, 1.9, 2.2, 5.15, 13.45, 5.85, 2.5, 5.3, 19.1]


def all_hotspot_and_acceptable_mutations_combinations_generator(all_mutations_summary_file_path):
    """
    Generator that yields all possible combinations of the hotspot and acceptable mutations.
    """
    all_mutations_summary_df = pd.read_csv(all_mutations_summary_file_path, header=0, index_col=0)

    acceptable_positions_mut_names_map = get_acceptable_positions_mut_names_map(all_mutations_summary_df)
    hotspot_positions_mut_names_map = get_hotspot_positions_mut_names_map(all_mutations_summary_df)
    
    # Remove hotspot positions from acceptable positions
    for position in hotspot_positions_mut_names_map.keys():
        del acceptable_positions_mut_names_map[position]

    ### NOTE: This only explores combinations of mutations in one order, do we want to get mutation order variability for BuildModel ? Would need to additionally use itertools.combinations.
    all_hotspot_mutations_combinations = list(itertools.product(*hotspot_positions_mut_names_map.values())) # hotspot_positions_mut_names_map.values() returns a list of sublists, where each sublist are the mutations for a given position
    if len(all_hotspot_mutations_combinations) == 0:
        all_hotspot_mutations_combinations = [[]]
    random.shuffle(all_hotspot_mutations_combinations)

    all_acceptable_mutations_combinations = list(itertools.product(*acceptable_positions_mut_names_map.values()))
    if len(all_acceptable_mutations_combinations) == 0:
        all_acceptable_mutations_combinations = [[]]
    random.shuffle(all_acceptable_mutations_combinations)

    for hotspot_mutations_list in all_hotspot_mutations_combinations:
        for acceptable_mutations_list in all_acceptable_mutations_combinations:
            mutations_list = hotspot_mutations_list + acceptable_mutations_list

            yield mutations_list
    
    return

def run_foldx_commands(mutant_PDB_file_path, wildtype_PDB_file_path, antibody_chains, antigen_chains, output_dir, GLOBALS):    
    run_foldx_AnalyseComplex(
        GLOBALS.foldx_dir,
        mutant_PDB_file_path.parent, mutant_PDB_file_path.stem,
        antibody_chains, antigen_chains,
        GLOBALS.vdwDesign,
        GLOBALS.print_stdout,
        output_dir, output_file_tag=None
    )
    run_foldx_AnalyseComplex(
        GLOBALS.foldx_dir,
        wildtype_PDB_file_path.parent, wildtype_PDB_file_path.stem,
        antibody_chains, antigen_chains,
        GLOBALS.vdwDesign,
        GLOBALS.print_stdout,
        output_dir, output_file_tag=None
    )

    run_foldx_Stability(
        GLOBALS.foldx_dir,
        mutant_PDB_file_path.parent, mutant_PDB_file_path.stem,
        GLOBALS.vdwDesign,
        GLOBALS.print_stdout,
        output_dir, output_file_tag=None
    )
    run_foldx_Stability(
        GLOBALS.foldx_dir,
        wildtype_PDB_file_path.parent, wildtype_PDB_file_path.stem,
        GLOBALS.vdwDesign,
        GLOBALS.print_stdout,
        output_dir, output_file_tag=None
    )

    if GLOBALS.calculate_binding_dG_with_water:
        run_foldx_AnalyseComplex(
            GLOBALS.foldx_dir,
            mutant_PDB_file_path.parent, mutant_PDB_file_path.stem,
            antibody_chains, antigen_chains,
            GLOBALS.vdwDesign,
            GLOBALS.print_stdout,
            output_dir, output_file_tag='mutant_with_waters',
            with_predicted_waters=True
        )
        run_foldx_AnalyseComplex(
            GLOBALS.foldx_dir,
            wildtype_PDB_file_path.parent, wildtype_PDB_file_path.stem,
            antibody_chains, antigen_chains,
            GLOBALS.vdwDesign,
            GLOBALS.print_stdout,
            output_dir, output_file_tag='wildtype_with_waters',
            with_predicted_waters=True
        )

    return

def systematic_search(parallel_executor, backbone_PDB_files_paths, evolvex_working_dir, GLOBALS):
    antibody_chains, antigen_chains = GLOBALS.antibody_chains, GLOBALS.antigen_chains

    for PDB_file_path in backbone_PDB_files_paths:
        PDB_dir = evolvex_working_dir / PDB_file_path.stem
        PDB_name = PDB_dir.name

        search_output_dir = PDB_dir / 'search_results'; search_output_dir.mkdir(exist_ok = True) ### Should be False
        all_mutations_summary_file_path = PDB_dir / 'hotspot_mutants' / 'all_mutations_summary.csv'

        futures = []
        for nth_model, mutations_list in enumerate(all_hotspot_and_acceptable_mutations_combinations_generator(all_mutations_summary_file_path)):
            foldx_Alanine_mutant_PDB_file_path = PDB_dir / f'{PDB_name}_1_Alanine_mutant.pdb'
            output_dir = search_output_dir / str(nth_model)

            future = parallel_executor.submit(
                create_model, foldx_Alanine_mutant_PDB_file_path, True, mutations_list, output_dir, GLOBALS,
            )
            futures.append(future)
        
        wait_and_remove(parallel_executor, futures)
    
    return


def get_random_mut_name(full_residue_IDs_list, allowed_AA_mutations_per_position_map):
    random_residue_ID = random.choice(full_residue_IDs_list)
    wildtype_AA, residue_ID, position = random_residue_ID[0], random_residue_ID[1:], random_residue_ID[2:]

    allowed_mutations = allowed_AA_mutations_per_position_map[position]
    if len(allowed_mutations) == 1:
        # Can't index a set, and there is only one mutation, so this works
        for mutant_AA in allowed_mutations:
            pass
    else:
        mutant_AA = random.choices(paratope_AA, paratope_AA_weights)[0]
        while not mutant_AA in allowed_mutations:
            mutant_AA = random.choices(paratope_AA, paratope_AA_weights)[0]
    
    return f'{wildtype_AA}{residue_ID}{mutant_AA}'

def metropolis_criterion(energies):
    if any(energy < 0 for energy in energies):
        return True
    
    U = random.random()
    for energy in energies:
        P = math.exp(-(energy / 0.5919))
        if P >= U:
            return True
        
    return False

def keep_mutant_decision(model_dir, antibody_chains, antigen_chains, antibody_stability_dG_original_wildtype, iteration_fraction, generated_models_info, GLOBALS):
    # Calculate energies and scores
    wildtype_antibody_stability_dG = get_chain_group_stability_dG(indiv_file_path = model_dir / 'Indiv_energies_WT_model_1_AC.fxout', chain_group_name = antibody_chains)
    mutant_antibody_stability_dG = get_chain_group_stability_dG(indiv_file_path = model_dir / 'Indiv_energies_model_1_AC.fxout', chain_group_name = antibody_chains)
    antibody_stability_ddG = round(mutant_antibody_stability_dG - wildtype_antibody_stability_dG, NDIGIS_ROUNDING)

    wildtype_complex_stability_dG = get_complex_stability_dG(st_file_path = model_dir / 'WT_model_1_0_ST.fxout')
    mutant_complex_stability_dG = get_complex_stability_dG(st_file_path = model_dir / 'model_1_0_ST.fxout')

    wildtype_binding_dG = get_binding_dG(interaction_file_path = model_dir / f'Interaction_WT_model_1_AC.fxout')
    mutant_binding_dG = get_binding_dG(interaction_file_path = model_dir / f'Interaction_model_1_AC.fxout')
    binding_ddG = round(mutant_binding_dG - wildtype_binding_dG, NDIGIS_ROUNDING)
    
    if GLOBALS.calculate_binding_dG_with_water:
        wildtype_binding_dG_with_waters = get_binding_dG(interaction_file_path = model_dir / f'Interaction_wildtype_with_waters_AC.fxout')
        mutant_binding_dG_with_waters = get_binding_dG(interaction_file_path = model_dir / f'Interaction_mutant_with_waters_AC.fxout')
        binding_ddG_with_waters = round(mutant_binding_dG_with_waters - wildtype_binding_dG_with_waters, NDIGIS_ROUNDING)

    wildtype_antibody_intraclash_score = get_chain_group_intraclash_score(interaction_file_path = model_dir / f'Interaction_WT_model_1_AC.fxout', chain_group_name = antibody_chains)
    mutant_antibody_intraclash_score = get_chain_group_intraclash_score(interaction_file_path = model_dir / f'Interaction_model_1_AC.fxout', chain_group_name = antibody_chains)
    antibody_delta_intraclash_score = round(mutant_antibody_intraclash_score - wildtype_antibody_intraclash_score, NDIGIS_ROUNDING)

    # All columns from interactions file, starting from Backbone Hbond
    wildtype_other_info_map = get_all_other_interaction_file_info(interaction_file_path = model_dir / f'Interaction_WT_model_1_AC.fxout')
    mutant_other_info_map = get_all_other_interaction_file_info(interaction_file_path = model_dir / f'Interaction_model_1_AC.fxout')

    
    # Make decision
    if antibody_stability_ddG > 0.5 or mutant_antibody_stability_dG > (antibody_stability_dG_original_wildtype + 2):
        keep_mutant =  False
    
    elif antibody_delta_intraclash_score > 0.5 or mutant_antibody_intraclash_score > 10:
        keep_mutant =  False
    
    else:
        energies = (binding_ddG, binding_ddG_with_waters) if GLOBALS.calculate_binding_dG_with_water else (binding_ddG,)
        keep_mutant = metropolis_criterion(energies)

    # Log info of the selected model. This bit of code is uggly, but couldn't find a better way.
    generated_models_info['antibody_stability_dG'].append(mutant_antibody_stability_dG if keep_mutant else wildtype_antibody_stability_dG)
    generated_models_info['complex_stability_dG'].append(mutant_complex_stability_dG if keep_mutant else wildtype_complex_stability_dG)
    generated_models_info['binding_dG'].append(mutant_binding_dG if keep_mutant else wildtype_binding_dG)
    if GLOBALS.calculate_binding_dG_with_water:
        generated_models_info['binding_dG_with_waters'].append(mutant_binding_dG_with_waters if keep_mutant else wildtype_binding_dG_with_waters)
    generated_models_info['antibody_intraclash_score'].append(mutant_antibody_intraclash_score if keep_mutant else wildtype_antibody_intraclash_score)
    
    other_info_map = mutant_other_info_map if keep_mutant else wildtype_other_info_map
    for key, value in other_info_map.items():
        key = key.replace(' ', '_') # Backbone Hbond => Backbone_Hbond
        generated_models_info[key].append(value)

    #

    return keep_mutant

def update_full_residue_IDs_list(model, mut_names_list):
    for mut_name in mut_names_list:
        wildtype_AA, residue_ID, mutant_AA = mut_name[0], mut_name[1:-1], mut_name[-1]
        old_residue_ID = f'{wildtype_AA}{residue_ID}'
        new_residue_ID = f'{mutant_AA}{residue_ID}'

        model.full_residue_IDs_list.remove(old_residue_ID)
        model.full_residue_IDs_list.append(new_residue_ID)

    return

def make_MC_steps(model, n_MC_steps, nth_loop, iteration_fraction, model_PDB_files_dir, GLOBALS):
    backbone_PDB_file_name = model.backbone_PDB_file_name
    antibody_stability_dG_original_wildtype = model.antibody_stability_dG_original_wildtype
    antibody_seq_map_original_wildtype = model.antibody_seq_map_original_wildtype
    allowed_AA_mutations_per_position_map = model.allowed_AA_mutations_per_position_map
    antibody_chains, antigen_chains = GLOBALS.antibody_chains, GLOBALS.antigen_chains

    generated_models_info = defaultdict(list)
    generated_models_info['backbone_PDB_file_name'] = [backbone_PDB_file_name] * n_MC_steps
    generated_models_info['nth_model'] = [model.model_dir.name] * n_MC_steps
    generated_models_info['step'] = ['MC'] * n_MC_steps

    current_nth_MC_iteration = nth_loop * (n_MC_steps + 1) # +1 to account for the recombination steps
    for i in range(n_MC_steps):
        nth_iteration = current_nth_MC_iteration + i + 1
        generated_models_info['nth_iteration'].append(nth_iteration)

        full_residue_IDs_list = model.full_residue_IDs_list
        mut_name = get_random_mut_name(full_residue_IDs_list, allowed_AA_mutations_per_position_map)

        # Create the mutant with BuildModel, which will be called "model_1"
        model_dir = model.model_dir
        create_model(
            input_PDB_file_path = model_dir / 'model.pdb', copy_PDB_file_to_output_dir = False, mutations_list = [mut_name], output_dir = model_dir, GLOBALS = GLOBALS,
        )

        run_foldx_commands(
            mutant_PDB_file_path = model_dir / 'model_1.pdb', wildtype_PDB_file_path = model_dir / 'WT_model_1.pdb',
            antibody_chains = antibody_chains, antigen_chains = antigen_chains, 
            output_dir = model_dir, GLOBALS = GLOBALS
        )

        keep_mutant = keep_mutant_decision(model_dir, antibody_chains, antigen_chains, antibody_stability_dG_original_wildtype, iteration_fraction, generated_models_info, GLOBALS)
        if keep_mutant:
            clean_up_model_dir(model_dir, PDB_file_name_to_keep_as_model = 'model_1.pdb')
            update_full_residue_IDs_list(model, mut_names_list = [mut_name])
        
        else:
            clean_up_model_dir(model_dir, PDB_file_name_to_keep_as_model = 'model.pdb')
        
        save_compressed_PDB_file(
            PDB_file_path = model_dir / 'model.pdb', 
            output_name = f'{backbone_PDB_file_name}_{model.model_dir.name}_{nth_iteration}.pdb',
            output_dir = model_PDB_files_dir
        )

        generated_models_info['residue_IDs'].append(';'.join(full_residue_IDs_list))
        generated_models_info['from_mut_name'].append(mut_name)
        generated_models_info['mutation_accepted'].append(keep_mutant)

    return model, generated_models_info

def write_generated_models_info(generated_models_info, generated_models_info_file_handle):
    if generated_models_info_file_handle.tell() == 0: # File is empty
        headers = ','.join(generated_models_info.keys()) + '\n'
        generated_models_info_file_handle.write(headers)

    # Equivalent to lines = pd.DataFrame(generated_models_info).to_csv(index=False, header=False), but 10x faster
    lines = ''.join(
        f"{','.join(map(str, row_of_values))}\n" # row CSV format
        for row_of_values in zip(*generated_models_info.values()) # if generated_models_info = {'A':[1,2,3], 'B':[4,5,6]}, zip(*generated_models_info.values()) yields (1,4), (2,5) and (3,6)
    )
    
    generated_models_info_file_handle.writelines(lines)
    return

def random_model_pairing_generator(models_population):
    backbone_PDB_grouped_models_map = defaultdict(list)
    for model in models_population:
        backbone_PDB_file_name = model.backbone_PDB_file_name
        backbone_PDB_grouped_models_map[backbone_PDB_file_name].append(model)

    for _, backbone_models_list in backbone_PDB_grouped_models_map.items():
        random.shuffle(backbone_models_list)

        for i in range(0, len(backbone_models_list), 2):
            yield (backbone_models_list[i], backbone_models_list[i+1])

    return

def get_recombination_mut_names(model_1, model_2):
    residue_ID_to_AA_map_1 = {full_residue_ID[1:]:full_residue_ID[0] for full_residue_ID in model_1.full_residue_IDs_list} # e.g {'H52':'K', 'H56':'Y', ...}
    residue_ID_to_AA_map_2 = {full_residue_ID[1:]:full_residue_ID[0] for full_residue_ID in model_2.full_residue_IDs_list}

    shared_residue_IDs = set.intersection(set(residue_ID_to_AA_map_1.keys()), set(residue_ID_to_AA_map_2.keys()))
    sorted_shared_residue_IDs = sorted(shared_residue_IDs, key = lambda residue_ID:(residue_ID[0], int(residue_ID[1:]))) 
    if len(sorted_shared_residue_IDs) < 2:
        raise ValueError(f'Cannot perform recombination between {model_1.model_dir} and {model_2.model_dir} because {shared_residue_IDs = }, and at least 2 shared residue IDs are needed.')

    recombination_location = random.randint(1, len(sorted_shared_residue_IDs) - 1)
    mut_names_1, mut_names_2 = [], []
    for residue_ID in sorted_shared_residue_IDs[:recombination_location]:
        AA_1 = residue_ID_to_AA_map_1[residue_ID]
        AA_2 = residue_ID_to_AA_map_2[residue_ID]

        mut_names_1.append(f'{AA_1}{residue_ID}{AA_2}') # e.g 'KH52R'
        mut_names_2.append(f'{AA_2}{residue_ID}{AA_1}')

    return (mut_names_1, mut_names_2)

def make_recombination_step(model_1, model_2, nth_iteration, iteration_fraction, model_PDB_files_dir, GLOBALS):
    # These variables are the same for both model 1 and 2
    backbone_PDB_file_name = model_1.backbone_PDB_file_name
    antibody_stability_dG_original_wildtype = model_1.antibody_stability_dG_original_wildtype
    antibody_seq_map_original_wildtype = model_1.antibody_seq_map_original_wildtype
    antibody_chains, antigen_chains = GLOBALS.antibody_chains, GLOBALS.antigen_chains

    generated_models_info = defaultdict(list)
    generated_models_info['backbone_PDB_file_name'] = [backbone_PDB_file_name] * 2
    generated_models_info['nth_model'] = [model_1.model_dir.name, model_2.model_dir.name]
    generated_models_info['step'] = ['recombination'] * 2
    generated_models_info['nth_iteration'] = [nth_iteration] * 2

    mut_names_1, mut_names_2 = get_recombination_mut_names(model_1, model_2)
    for mut_names, model in [(mut_names_1, model_1), (mut_names_2, model_2)]:
        model_dir = model.model_dir
        full_residue_IDs_list = model.full_residue_IDs_list

        create_model(
            input_PDB_file_path = model_dir / 'model.pdb', copy_PDB_file_to_output_dir = False, mutations_list = mut_names, output_dir = model_dir, GLOBALS = GLOBALS
        )

        run_foldx_commands(
            mutant_PDB_file_path = model_dir / 'model_1.pdb', wildtype_PDB_file_path = model_dir / 'WT_model_1.pdb', 
            antibody_chains = antibody_chains, antigen_chains = antigen_chains, 
            output_dir = model_dir, GLOBALS = GLOBALS
        )
        
        keep_mutant = keep_mutant_decision(model_dir, antibody_chains, antigen_chains, antibody_stability_dG_original_wildtype, iteration_fraction, generated_models_info, GLOBALS)
        if keep_mutant:
            clean_up_model_dir(model_dir, PDB_file_name_to_keep_as_model = 'model_1.pdb')
            update_full_residue_IDs_list(model, mut_names_list = mut_names)
        
        else:
            clean_up_model_dir(model_dir, PDB_file_name_to_keep_as_model = 'model.pdb')

        save_compressed_PDB_file(
            PDB_file_path = model_dir / 'model.pdb', 
            output_name = f'{backbone_PDB_file_name}_{model.model_dir.name}_{nth_iteration}.pdb',
            output_dir = model_PDB_files_dir
        )

        generated_models_info['residue_IDs'].append(';'.join(full_residue_IDs_list))
        generated_models_info['from_mut_name'].append(';'.join(mut_names))
        generated_models_info['mutation_accepted'].append(keep_mutant)
    
    return model_1, model_2, generated_models_info

def GA_search(parallel_executor, initial_models_population, generated_models_info_file_path, model_PDB_files_dir, GLOBALS):
    recombine_every_nth_iteration = GLOBALS.recombine_every_nth_iteration
    max_iterations = GLOBALS.max_iterations

    # Each loop = n MC steps + 1 recombination step
    n_MC_and_recombination_loops = max_iterations // recombine_every_nth_iteration
    n_MC_steps_per_loop = recombine_every_nth_iteration - 1

    generated_models_info_file_handle = open(generated_models_info_file_path, 'a')

    models_population = initial_models_population
    for nth_loop in range(n_MC_and_recombination_loops):
        iteration_fraction = nth_loop / n_MC_and_recombination_loops

        # MC steps
        futures = []
        for model in models_population:
            future = parallel_executor.submit(make_MC_steps, model, n_MC_steps_per_loop, nth_loop, iteration_fraction, model_PDB_files_dir, GLOBALS)
            futures.append(future)

        models_population = []
        for _, (model, generated_models_info) in as_completed(futures, with_results=True):
            models_population.append(model)
            write_generated_models_info(generated_models_info, generated_models_info_file_handle)
        parallel_executor.cancel(futures)
        
        # Recombination step. The models are recombined with a model from the same backbone.
        nth_iteration = (nth_loop + 1) * recombine_every_nth_iteration
        futures = []
        for model_1, model_2 in random_model_pairing_generator(models_population):
            future = parallel_executor.submit(make_recombination_step, model_1, model_2, nth_iteration, iteration_fraction, model_PDB_files_dir, GLOBALS)
            futures.append(future)

        models_population = []
        for _, (model_1, model_2, generated_models_info) in as_completed(futures, with_results=True):
            models_population.append(model_1)
            models_population.append(model_2)
            write_generated_models_info(generated_models_info, generated_models_info_file_handle)
        parallel_executor.cancel(futures)

    generated_models_info_file_handle.close()
    return