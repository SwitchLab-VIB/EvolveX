import warnings

from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1

from Bio import BiopythonParserWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter('ignore', PDBConstructionWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)

def get_residue_ID_to_residue_name_map(PDB_file_path):
    """
    Returns a dict where keys are residue IDs (e.g A52, for residue 52 in chain A) and keys are single letter residues (e.g G for Glycine)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(id='', file=PDB_file_path)

    residue_ID_to_residue_name_map = {}
    for residue in structure.get_residues():
        _, _, chain, (_, number, _) = residue.full_id
        residue_ID = f'{chain}{number}'
        residue_name = protein_letters_3to1[residue.resname.title()]

        residue_ID_to_residue_name_map[residue_ID] = residue_name

    return residue_ID_to_residue_name_map

def get_chain_to_sequence_map(PDB_file_path, chain_subset):
    chain_to_sequence_map = {
        record.annotations['chain']:str(record.seq).replace('X', '') # When extracting the sequences from ATOMS lines, SeqIO.parse introduces an 'X' when there are consecutive residues with a number difference > 1 (e.g: A1,Y2,P5 => AYXXP), but we don't want that.
        for record in SeqIO.parse(PDB_file_path, format = 'pdb-atom')
        if record.annotations['chain'] in chain_subset
    }
    return chain_to_sequence_map