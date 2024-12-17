
import tarfile

NDIGIS_ROUNDING = 4

def save_compressed_PDB_file(PDB_file_path, output_name, output_dir):
    with tarfile.open(output_dir / f'{output_name}.tar.gz', 'w:gz') as tar_file_handle:
        tar_file_handle.add(PDB_file_path, arcname=output_name)
    return