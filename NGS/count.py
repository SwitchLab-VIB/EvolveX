import os
import sys
import subprocess
from pathlib import Path

# ---------------- Setup virtual environment ---------------- #
venv_dir = Path(__file__).parent / "venv"
venv_python = venv_dir / "bin" / "python"

if not venv_dir.exists():
    print("🔧 Creating virtual environment...")
    subprocess.run([sys.executable, "-m", "venv", str(venv_dir)])
    print("✅ Virtual environment created.")

# Install dependencies inside venv if missing
def ensure_packages():
    try:
        import pandas  # noqa
        import Bio     # noqa
    except ImportError:
        print("📦 Installing required packages (pandas, biopython)...")
        subprocess.run([str(venv_python), "-m", "pip", "install", "--upgrade", "pip"])
        subprocess.run([str(venv_python), "-m", "pip", "install", "pandas", "biopython"])
        print("✅ Packages installed.")

ensure_packages()

# Re-run inside venv if not already
if sys.executable != str(venv_python):
    print("🔁 Re-running inside virtual environment...")
    os.execv(str(venv_python), [str(venv_python), __file__])

# ---------------- Actual imports ---------------- #
import pandas as pd
from Bio.Seq import Seq
import time
import gc

# ---------------- Processing logic ---------------- #
output_summary_file = "all_processing_summary.csv"

# Remove existing summary file
if os.path.exists(output_summary_file):
    os.remove(output_summary_file)

for i in range(63, 86):
    file_name = f"n{i}.tabular"
    if not os.path.exists(file_name):
        print(f"⚠️ Skipping {file_name} (not found)")
        continue

    start_time = time.time()
    print(f"\n▶️ Processing {file_name} ...")

    # Read file
    df = pd.read_csv(file_name, sep='\t', header=None, names=['sequence'])

    # Extract UMI and sequence info
    df['UMI-Fw'] = df['sequence'].str[:15]
    df['UMI-Rv'] = df['sequence'].str[-15:]
    df['UMI pair'] = df['UMI-Fw'] + df['UMI-Rv']

    reads_before = len(df)
    unique_umi_pair_before = df['UMI pair'].nunique()
    unique_umi_fw_before = df['UMI-Fw'].nunique()
    unique_umi_rv_before = df['UMI-Rv'].nunique()

    df.drop_duplicates(inplace=True)
    after_duplicates = len(df)

    df['extracted_sequence'] = df['sequence'].str[34:229]
    df['amino_acid_sequence'] = df['extracted_sequence'].apply(lambda x: str(Seq(x).translate()))
    df['amino_acid_sequence'] = df['amino_acid_sequence'].str.replace(',', '', regex=False)

    df['CDR2'] = df['amino_acid_sequence'].str[:10]
    df['CDR3'] = df['amino_acid_sequence'].str[-16:]

    df = df[~df['amino_acid_sequence'].str.contains(r'\*')]
    after_stop_codons = len(df)

    # Save counts
    pd.DataFrame(df['CDR2'].value_counts()).reset_index() \
        .rename(columns={'index': 'CDR2', 'CDR2':'Count_CDR2'}) \
        .to_csv(f'n{i}_cdr2_counts.csv', index=False)
    pd.DataFrame(df['CDR3'].value_counts()).reset_index() \
        .rename(columns={'index': 'CDR3', 'CDR3':'Count_CDR3'}) \
        .to_csv(f'n{i}_cdr3_counts.csv', index=False)

    # Prepare summary
    summary_data = {
        "File": f"n{i}",
        "Number of reads before": reads_before,
        "Number of unique UMI pair before": unique_umi_pair_before,
        "Number of unique UMI-Fw before": unique_umi_fw_before,
        "Number of unique UMI-Rv before": unique_umi_rv_before,
        "Number of duplicates removed": reads_before - after_duplicates,
        "Number of reads left after removing duplicates": after_duplicates,
        "Number of stop codons removed": after_duplicates - after_stop_codons,
        "Number of reads left after removing stop codons": after_stop_codons,
        "Number of unique UMI pair after processing": df['UMI pair'].nunique()
    }

    # Append summary to CSV
    summary_df = pd.DataFrame([summary_data])
    summary_df.to_csv(output_summary_file, mode='a', index=False, header=not os.path.exists(output_summary_file))

    # Free memory
    del df
    gc.collect()

    print(f"✅ Finished {file_name} in {time.time() - start_time:.2f} seconds")

print(f"\n📊 All summaries merged into {output_summary_file}")
