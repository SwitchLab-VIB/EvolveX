import pandas as pd

from Bio.Data.IUPACData import protein_letters

from dask.distributed import as_completed

from evolvex.dask_parallel import setup_dask_parallel_executor
from evolvex.mutate_interface import generate_Alanine_mutant, mutate_antibody_hotspot_position, generate_mutations_summary_file
from evolvex.search_algorithms import GA_search, systematic_search
from evolvex.model_generation import generate_initial_models
from evolvex.dask_parallel import wait_and_remove
from evolvex.command_line_interface import command_line_interface

def main():
    """
    """
    GLOBALS = command_line_interface()
    
    evolvex_working_dir = GLOBALS.working_dir / 'EvolveX'; evolvex_working_dir.mkdir(exist_ok = True) ### Should be false

    backbone_PDB_files_paths = list(GLOBALS.Backbones_dir.glob('*.pdb'))
    all_PDBs_positions_to_explore_df = pd.read_csv(GLOBALS.PositionsToExplore_file_path, header=0, sep='\t')
    
    parallel_executor = setup_dask_parallel_executor(GLOBALS)
    

    # Generate Alanine mutants of all PDB backbones
    futures_1 = []
    for PDB_file_path in backbone_PDB_files_paths:
        PDB_positions_to_explore_df = all_PDBs_positions_to_explore_df[all_PDBs_positions_to_explore_df['Pdb'] == PDB_file_path.stem]

        future_1 = parallel_executor.submit(
            generate_Alanine_mutant, PDB_file_path, PDB_positions_to_explore_df, evolvex_working_dir, GLOBALS
        )
        futures_1.append(future_1)

    # Explore all possible mutations for positions marked as 'AUTO'
    print('Exploring mutations at each position...', flush=True)
    residues_to_explore = set(protein_letters) - set(GLOBALS.residues_to_ignore)
    futures_2 = []
    for future_1, result in as_completed(futures_1, with_results=True):
        foldx_Alanine_mutant_PDB_file_path, AUTO_Ala_positions_full_residue_IDs, output_dir = result
        future_1.release()
        for full_residue_ID in AUTO_Ala_positions_full_residue_IDs:
            hotspot_mutants_dir = output_dir / 'hotspot_mutants' / full_residue_ID
            hotspot_mutants_dir.mkdir(parents=True, exist_ok=True)
            for mutant_residue in residues_to_explore:
                future_2 = parallel_executor.submit(
                    mutate_antibody_hotspot_position, foldx_Alanine_mutant_PDB_file_path, full_residue_ID, mutant_residue, hotspot_mutants_dir, GLOBALS
                )
                futures_2.append(future_2)

    wait_and_remove(parallel_executor, futures_2)

    # Generate a summary file for each PDB backbone, which includes the ddG_binding, ddG_stability_complex and ddG_stability_antibody for all possible 
    # mutations at each position. For positions with pre-selected mutations, the fields are artificially set to -100.
    futures_3 = []
    for PDB_file_path in backbone_PDB_files_paths:
        PDB_dir = evolvex_working_dir / PDB_file_path.stem
        PDB_positions_to_explore_df = all_PDBs_positions_to_explore_df[all_PDBs_positions_to_explore_df['Pdb'] == PDB_file_path.stem]

        future_3 = parallel_executor.submit(
            generate_mutations_summary_file, PDB_dir, PDB_positions_to_explore_df, GLOBALS
        )
        futures_3.append(future_3)

    wait_and_remove(parallel_executor, futures_3)

    # Run search algorithm
    if GLOBALS.search_algorithm == 'systematic':
        systematic_search(parallel_executor, backbone_PDB_files_paths, evolvex_working_dir, GLOBALS)
    
    else:
        generated_models_info_file_path = GLOBALS.working_dir / 'generated_models_info.csv'
        model_PDB_files_dir = GLOBALS.working_dir / 'model_PDB_files'; model_PDB_files_dir.mkdir(exist_ok = True)

        print('Generating initial models...', flush=True)
        initial_models_population = generate_initial_models(parallel_executor, evolvex_working_dir, backbone_PDB_files_paths, GLOBALS)

        print('Running search...', flush=True)
        GA_search(parallel_executor, initial_models_population, generated_models_info_file_path, model_PDB_files_dir, GLOBALS)


    print('Finished.', flush=True)
    parallel_executor.close()
    return

if __name__ == '__main__':
    main()