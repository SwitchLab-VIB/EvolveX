import shutil

import pandas as pd

from evolvex.foldx_commands import (
    create_individual_list_foldx_mutations_file, get_alanine_mutant,
    run_foldx_BuildModel, run_foldx_AnalyseComplex, run_foldx_Stability,
    get_binding_ddG, get_complex_stability_ddG, get_chain_group_stability_ddG
)

def generate_Alanine_mutant(PDB_file_path, PDB_positions_to_explore_df, evolvex_working_dir, GLOBALS):
    PDB_file_name = PDB_file_path.stem
    output_dir = evolvex_working_dir / PDB_file_name; output_dir.mkdir(exist_ok=True)

    # PDB_file_path corresponds to the original file in the Backbones folder, which we do not want to modify. The PDB_file_path_copy is 
    # located in the corresponding subfolder named after the PDB file name in evolvex_working_dir
    PDB_file_path_copy = shutil.copy(
        src = PDB_file_path, 
        dst = output_dir / f'{PDB_file_name}.pdb'
    )

    antibody_chains, antigen_chains = GLOBALS.antibody_chains, GLOBALS.antigen_chains

    # PositionsToExplore possibilities:
    # 1) AA_Allowed=AUTO & MakeAla=Y --> Mutate to Alanine and find which mutations are worth exploring automatically 
    # 2) AA_Allowed=string of AA & MakeAla=Y --> Mutate to Alanine and only allow the specified mutations
    # 3) AA_Allowed=string of AA & MakeAla=N --> Do not mutate to Alanine (i.e keep the wildtype residue) and only allow the specified mutations
    PDB_positions_to_explore_df['full_residue_ID'] = PDB_positions_to_explore_df.apply(lambda row:f'{row.Res1}{row.Chain}{row.number}', axis=1)
    positions_to_Ala_df = PDB_positions_to_explore_df[PDB_positions_to_explore_df['MakeAla'] == 'Y']
    if positions_to_Ala_df.empty:
        raise ValueError('There must be at least 1 position to mutate to Alanine.')

    positions_to_Ala_full_residue_IDs = positions_to_Ala_df.full_residue_ID.values
    Alanine_mutant = get_alanine_mutant(positions_to_Ala_full_residue_IDs)
    individual_list_foldx_mutations_file_path = create_individual_list_foldx_mutations_file(mutant = Alanine_mutant, output_dir = output_dir)

    foldx_Alanine_mutant_PDB_file_path, foldx_wildtype_PDB_file_path = run_foldx_BuildModel(
        foldx_dir=GLOBALS.foldx_dir, 
        PDB_file_dir=PDB_file_path_copy.parent, PDB_file_name=PDB_file_name, 
        individual_list_foldx_mutations_file_path=individual_list_foldx_mutations_file_path, 
        move_neighbors_flag=False, # At this stage we do not need an optimized structure, it will get optimized latter
        vdwDesign=GLOBALS.vdwDesign,
        print_stdout=GLOBALS.print_stdout,
        output_dir=output_dir, output_file_tag='Alanine_mutant', PDB_file_tag='Alanine_mutant'
    )

    # The stability of the original wildtype antibody is needed for the MC and GA searches
    run_foldx_AnalyseComplex(
        GLOBALS.foldx_dir,
        PDB_file_path_copy.parent, PDB_file_path_copy.stem,
        antibody_chains, antigen_chains,
        GLOBALS.vdwDesign,
        GLOBALS.print_stdout,
        output_dir, output_file_tag='original_wildtype'
    )

    AUTO_Ala_positions_full_residue_IDs = positions_to_Ala_df[positions_to_Ala_df['AA_Allowed'] == 'AUTO'].full_residue_ID.values
    return foldx_Alanine_mutant_PDB_file_path, AUTO_Ala_positions_full_residue_IDs, output_dir 

def mutate_antibody_hotspot_position(foldx_Alanine_mutant_PDB_file_path, full_residue_ID, mutant_residue, hotspot_mutants_dir, GLOBALS):
    antibody_chains, antigen_chains = GLOBALS.antibody_chains, GLOBALS.antigen_chains
    hotspot_mutant = f'A{full_residue_ID[1:]}{mutant_residue}' # e.g: AH25K

    hotspot_mutant_output_dir = hotspot_mutants_dir / hotspot_mutant; hotspot_mutant_output_dir.mkdir(exist_ok=True)
    shutil.copy(
        src = foldx_Alanine_mutant_PDB_file_path, 
        dst = hotspot_mutant_output_dir
    ) 
    
    individual_list_foldx_mutations_file_path = create_individual_list_foldx_mutations_file(mutant = hotspot_mutant, output_dir = hotspot_mutant_output_dir)
    foldx_mutant_PDB_file_path, foldx_wildtype_PDB_file_path = run_foldx_BuildModel(
        foldx_dir=GLOBALS.foldx_dir,
        PDB_file_dir=hotspot_mutant_output_dir, PDB_file_name=foldx_Alanine_mutant_PDB_file_path.stem,
        individual_list_foldx_mutations_file_path=individual_list_foldx_mutations_file_path,
        move_neighbors_flag=True, # At this stage we want to optimize the generated structures
        vdwDesign=GLOBALS.vdwDesign,
        print_stdout=GLOBALS.print_stdout,
        output_dir=hotspot_mutant_output_dir, output_file_tag='hotspot_mutant', PDB_file_tag='hotspot_mutant'
    )

    # Mutant
    run_foldx_AnalyseComplex(
        GLOBALS.foldx_dir,
        hotspot_mutant_output_dir, foldx_mutant_PDB_file_path.stem,
        antibody_chains, antigen_chains,
        GLOBALS.vdwDesign,
        GLOBALS.print_stdout,
        hotspot_mutant_output_dir, output_file_tag='hotspot_mutant'
    )
    run_foldx_Stability(
        GLOBALS.foldx_dir,
        hotspot_mutant_output_dir, foldx_mutant_PDB_file_path.stem,
        GLOBALS.vdwDesign,
        GLOBALS.print_stdout,
        hotspot_mutant_output_dir, output_file_tag='hotspot_mutant'
    )

    # Wildtype
    run_foldx_AnalyseComplex(
        GLOBALS.foldx_dir,
        hotspot_mutant_output_dir, foldx_wildtype_PDB_file_path.stem,
        antibody_chains, antigen_chains,
        GLOBALS.vdwDesign,
        GLOBALS.print_stdout,
        hotspot_mutant_output_dir, output_file_tag='hotspot_wildtype'
    )
    run_foldx_Stability(
        GLOBALS.foldx_dir,
        foldx_wildtype_PDB_file_path.parent, foldx_wildtype_PDB_file_path.stem,
        GLOBALS.vdwDesign,
        GLOBALS.print_stdout,
        hotspot_mutant_output_dir, output_file_tag='hotspot_wildtype'
    )

    return

def generate_mutations_summary_file(PDB_dir, PDB_positions_to_explore_df, GLOBALS):
    """
    """
    PDB_name = PDB_dir.name

    # Create a summary_df per position, and then combine them all into one final df
    all_mutations_summary_df = []
    for _, row in PDB_positions_to_explore_df.iterrows():
        full_original_residue_ID = f'{row.Res1}{row.Chain}{row.number}'

        output_dir = PDB_dir / 'hotspot_mutants' / full_original_residue_ID

        binding_ddG_position, complex_stability_ddG_position, antibody_stability_ddG_position, index = [], [], [], []
        if row.AA_Allowed != 'AUTO':
            assert len(row.AA_Allowed) > 0, f"Rows where AA_Allowed is not 'AUTO' must have at least one specified residue ({row = })"
            
            output_dir.mkdir(exist_ok=True)

            # Generate an 'artificial' summary df where all fields are set to -100 for each amino acid
            values = [-100] * len(row.AA_Allowed)
            binding_ddG_position.extend(values)
            complex_stability_ddG_position.extend(values)
            antibody_stability_ddG_position.extend(values)

            wildtype_AA = 'A' if row.MakeAla == 'Y' else full_original_residue_ID[0]
            index.extend((f'{wildtype_AA}{full_original_residue_ID[1:]}{mutant_AA}' for mutant_AA in row.AA_Allowed))

        else:
            for mut_name_dir in output_dir.iterdir(): # Iterate over each individual mutation folder that was generated (e.g AH27C, AH27D, ...) 
                if not mut_name_dir.is_dir():
                    continue
            
                binding_ddG = get_binding_ddG(
                    wildtype_interaction_file_path = mut_name_dir / 'Interaction_hotspot_wildtype_AC.fxout', # Wildtype = alanine mutant
                    mutant_interaction_file_path = mut_name_dir / 'Interaction_hotspot_mutant_AC.fxout'
                )
                complex_stability_ddG = get_complex_stability_ddG(
                    dif_file_path = mut_name_dir / f'Dif_hotspot_mutant_{PDB_name}_1_Alanine_mutant.fxout'
                )
                antibody_stability_ddG = get_chain_group_stability_ddG(
                    wildtype_indiv_file_path = mut_name_dir / 'Indiv_energies_hotspot_wildtype_AC.fxout',
                    mutant_indiv_file_path = mut_name_dir / 'Indiv_energies_hotspot_mutant_AC.fxout',
                    chain_group_name = GLOBALS.antibody_chains
                )

                binding_ddG_position.append(binding_ddG)
                complex_stability_ddG_position.append(complex_stability_ddG)
                antibody_stability_ddG_position.append(antibody_stability_ddG)
                
                index.append(mut_name_dir.name)


        summary_df = pd.DataFrame(
            data = {
                'binding_ddG':binding_ddG_position,
                'complex_stability_ddG':complex_stability_ddG_position,
                'antibody_stability_ddG':antibody_stability_ddG_position
            },
            index = index
        ).sort_index()
        summary_df['original_residue'] = full_original_residue_ID[0]
        summary_df['position'] = full_original_residue_ID[2:]

        all_mutations_summary_df.append(summary_df)

    if len(all_mutations_summary_df) == 0:
        print(f'Something went wrong when generating the mutations summary file for {PDB_name = }.')
        
    pd.concat(all_mutations_summary_df, axis=0).round(decimals = 6).to_csv(PDB_dir / 'hotspot_mutants' / 'all_mutations_summary.csv') 
    return