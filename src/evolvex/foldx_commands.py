import shutil
import subprocess
import shutil

from evolvex.utils import NDIGIS_ROUNDING

def get_alanine_mutant(full_residue_IDs):
    """
    Given the list of full_residue_IDs, returns a comma separated string of mutation names to Alanine, 
    following the FoldX naming convention (e.g 'LH4A,TH9A,KH17A').
    """
    mutation_names = (
        f'{full_residue_ID}A' # e.g KH52A
        for full_residue_ID in full_residue_IDs
    )
    return ','.join(mutation_names)

def create_individual_list_foldx_mutations_file(mutant, output_dir, output_file_name='individual_list_foldx_mutations_file.txt'):
    """
    Mutant should be a comma separated string of mutations that describe a mutant (e.g: "DH52A,LH80A,KH99A").

    By default, generates a file named "individual_list_foldx_mutations_file.txt" in output_dir.
    """
    output_file_path = output_dir / output_file_name
    with open(output_file_path, 'wt') as file_handle:
        file_handle.write(f'{mutant};\n')

    return output_file_path


def run_foldx_BuildModel(
        foldx_dir, PDB_file_dir, PDB_file_name, individual_list_foldx_mutations_file_path, move_neighbors_flag, vdwDesign, print_stdout, output_dir, output_file_tag, PDB_file_tag=None
    ):
    """
    Returns the full paths of the mutant and wildtype PDB file generated by BuildModel.
    """
    input_file = PDB_file_dir / f'{PDB_file_name}.pdb'
    assert input_file.exists(), f'{input_file} does not exist.'
    
    command = [
        str(foldx_dir / 'foldx'), 
        '--command', 'BuildModel',
        '--pdb-dir', str(PDB_file_dir),
        '--pdb', f'{PDB_file_name}.pdb',
        '--mutant-file', str(individual_list_foldx_mutations_file_path),
        '--pdbHydrogens', 'true',
        '--output-dir', str(output_dir),
        '--moveNeighbours', 'true' if move_neighbors_flag else 'false',
        '--vdwDesign', str(vdwDesign),
        '--screen', 'true' if print_stdout else 'false'
    ]
    if output_file_tag:
        command += ['--output-file', output_file_tag]

    subprocess.run(command, check=True, stdout=None if print_stdout else subprocess.DEVNULL)
    
    foldx_mutant_PDB_file_path   = PDB_file_dir / f'{PDB_file_name}_1.pdb'
    foldx_wildtype_PDB_file_path = PDB_file_dir / f'WT_{PDB_file_name}_1.pdb'

    # A new PDB of the wildtype is not generated by FoldX when move_neighbors_flag is False, as it simply corresponds to the input PDB file,
    # but to keep it consistent we create a copy of the file with the expected name.
    if move_neighbors_flag == False:
        shutil.copy(
            src = PDB_file_dir / f'{PDB_file_name}.pdb', 
            dst = foldx_wildtype_PDB_file_path
        )
    
    # As for the other files generated by BuildModel, gives the possibility to add a tag to the generated PDB files.
    if PDB_file_tag:
        foldx_mutant_PDB_file_path = foldx_mutant_PDB_file_path.rename(foldx_mutant_PDB_file_path.with_stem(f"{foldx_mutant_PDB_file_path.stem}_{PDB_file_tag}"))
        foldx_wildtype_PDB_file_path = foldx_wildtype_PDB_file_path.rename(foldx_wildtype_PDB_file_path.with_stem(f"{foldx_wildtype_PDB_file_path.stem}_{PDB_file_tag}"))

    return (foldx_mutant_PDB_file_path, foldx_wildtype_PDB_file_path)

def run_foldx_AnalyseComplex(foldx_dir, PDB_file_dir, PDB_file_name, antibody_chains, antigen_chains, vdwDesign, print_stdout, output_dir, output_file_tag, with_predicted_waters=False):
    input_file = PDB_file_dir / f'{PDB_file_name}.pdb'
    assert input_file.exists(), f'{input_file} does not exist.'

    command = [
        str(foldx_dir / 'foldx'), 
        '--command', 'AnalyseComplex',
        '--pdb-dir', str(PDB_file_dir),
        '--pdb', f'{PDB_file_name}.pdb',
        '--analyseComplexChains', f'{antibody_chains},{antigen_chains}',
        '--vdwDesign', str(vdwDesign),
        '--output-dir', str(output_dir),
        '--screen', 'true' if print_stdout else 'false'
    ]
    
    if with_predicted_waters:
        command += [
            '--water', '-PREDICT',
            '--ionStrength', '0.150',
        ]

    if output_file_tag:
        command += ['--output-file', output_file_tag]

    subprocess.run(command, check=True, stdout=None if print_stdout else subprocess.DEVNULL)
    return

def run_foldx_Stability(foldx_dir, PDB_file_dir, PDB_file_name, vdwDesign, print_stdout, output_dir, output_file_tag):
    input_file = PDB_file_dir / f'{PDB_file_name}.pdb'
    assert input_file.exists(), f'{input_file} does not exist.'

    command = [
        str(foldx_dir / 'foldx'), 
        '--command', 'Stability',
        '--pdb-dir', str(PDB_file_dir),
        '--pdb', f'{PDB_file_name}.pdb',
        '--vdwDesign', str(vdwDesign),
        '--output-dir', str(output_dir),
    ]
    if output_file_tag:
        command += ['--output-file', output_file_tag]

    subprocess.run(command, check=True, stdout=None if print_stdout else subprocess.DEVNULL)
    return


def get_binding_dG(interaction_file_path):
    assert interaction_file_path.name.startswith('Interaction_'), "To obtain binding dG, provide an 'Interaction_' file generated by AnalyseComplex."

    with open(interaction_file_path, 'rt') as file_handle:
        lines = file_handle.readlines()

    binding_dG = float( lines[9].split('\t')[5] )
    return round(binding_dG, NDIGIS_ROUNDING)

def get_binding_ddG(wildtype_interaction_file_path, mutant_interaction_file_path):
    assert wildtype_interaction_file_path.name.startswith('Interaction_') and mutant_interaction_file_path.name.startswith('Interaction_'), "To obtain binding ddG, provide two 'Interaction_' files generated by AnalyseComplex."

    wildtype_binding_dG = get_binding_dG(wildtype_interaction_file_path)
    mutant_binding_dG = get_binding_dG(mutant_interaction_file_path)

    binding_ddG = mutant_binding_dG - wildtype_binding_dG
    return round(binding_ddG, NDIGIS_ROUNDING)

def get_chain_group_stability_dG(indiv_file_path, chain_group_name):
    assert indiv_file_path.name.startswith('Indiv_'), "To obtain the stability dG of a chain group, provide an 'Indiv_' file generated by AnalyseComplex"
    
    with open(indiv_file_path, 'rt') as file_handle:
        lines = file_handle.readlines()

    for line in lines[9:]:
        line = line.split('\t')
        if line[1] == chain_group_name:
            chain_group_stability_dG = float( line[2] )
            return round(chain_group_stability_dG, NDIGIS_ROUNDING)
        
    raise ValueError(f'Could not find {chain_group_name=} in {indiv_file_path}')
    return

def get_chain_group_stability_ddG(wildtype_indiv_file_path, mutant_indiv_file_path, chain_group_name):
    assert wildtype_indiv_file_path.name.startswith('Indiv_') and mutant_indiv_file_path.name.startswith('Indiv_'), "To obtain the stability ddG of a chain group, provide two 'Indiv_' files generated by AnalyseComplex."

    wildtype_chain_group_stability_dG = get_chain_group_stability_dG(wildtype_indiv_file_path, chain_group_name)
    mutant_chain_group_stability_dG = get_chain_group_stability_dG(mutant_indiv_file_path, chain_group_name)

    chain_group_stability_ddG = mutant_chain_group_stability_dG - wildtype_chain_group_stability_dG
    return round(chain_group_stability_ddG, NDIGIS_ROUNDING)

def get_complex_stability_dG(st_file_path):
    # NOTE: The same information could be obtained from the 'Raw_' file generated by BuildModel, so we could skip running Stability completely
    assert st_file_path.name.endswith('_ST.fxout'), "To obtain the stability dG of a complex, provide a '_ST.fxout' file generated by Stability."

    with open(st_file_path, 'rt') as file_handle:
        lines = file_handle.readlines()

    complex_stability_dG = float( lines[0].split('\t')[1] )
    return round(complex_stability_dG, NDIGIS_ROUNDING)

def get_complex_stability_ddG(dif_file_path):
    assert dif_file_path.name.startswith('Dif_'), "To obtain the stability ddG of a complex, provide a 'Dif_' file generated by BuildModel."
    
    with open(dif_file_path, 'rt') as file_handle:
        lines = file_handle.readlines()

    complex_stability_ddG = float( lines[9].split('\t')[1] )
    return round(complex_stability_ddG, NDIGIS_ROUNDING)

def get_chain_group_intraclash_score(interaction_file_path, chain_group_name):
    assert interaction_file_path.name.startswith('Interaction_'), "To obtain intraclash scores of a chain group, provide an 'Interaction_' file generated by AnalyseComplex."

    with open(interaction_file_path, 'rt') as file_handle:
        lines = file_handle.readlines()
    
    line = lines[9].split('\t')
    intraclash_scores = {
        line[1]:float(line[3]), # e.g 'HL':2.363
        line[2]:float(line[4])
    }
    return round(intraclash_scores[chain_group_name], NDIGIS_ROUNDING)

def get_chain_group_delta_intraclash_score(wildtype_interaction_file_path, mutant_interaction_file_path, chain_group_name):
    assert wildtype_interaction_file_path.name.startswith('Interaction_') and mutant_interaction_file_path.name.startswith('Interaction_'), "To obtain a change in intraclash score of a chain group, provide two 'Interaction_' files generated by AnalyseComplex."
    
    wildtype_intraclash_score = get_chain_group_intraclash_score(wildtype_interaction_file_path, chain_group_name)
    mutant_intraclash_score = get_chain_group_intraclash_score(mutant_interaction_file_path, chain_group_name)
    
    chain_group_delta_intraclash_score = mutant_intraclash_score - wildtype_intraclash_score
    return round(chain_group_delta_intraclash_score, NDIGIS_ROUNDING)

def get_all_other_interaction_file_info(interaction_file_path):
    assert interaction_file_path.name.startswith('Interaction_'), "To obtain all the information to the right of 'Backbone Hbond' from an interaction file, provide an 'Interaction_' file generated by AnalyseComplex."

    with open(interaction_file_path, 'rt') as file_handle:
        lines = file_handle.readlines()

    column_names = lines[8].strip().split('\t')
    values = lines[9].strip().split('\t')

    return {column_names[i]:round(float(values[i]), NDIGIS_ROUNDING) for i in range(6, len(column_names))}
