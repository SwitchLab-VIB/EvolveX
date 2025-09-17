import os
import sys
import subprocess
from pathlib import Path

# Step 1: Check for virtual environment
venv_dir = Path(__file__).parent / "venv"
venv_python = venv_dir / "bin" / "python"

if not venv_dir.exists():
    print("🔧 Creating virtual environment...")
    subprocess.run([sys.executable, "-m", "venv", str(venv_dir)])
    print("✅ Virtual environment created.")

# Step 2: Check if dependencies are installed, if not install them
try:
    import pandas
    import Bio
except ImportError:
    print("📦 Installing required packages (pandas, biopython)...")
    subprocess.run([str(venv_python), "-m", "pip", "install", "--upgrade", "pip"])
    subprocess.run([str(venv_python), "-m", "pip", "install", "pandas", "biopython"])
    print("✅ Packages installed.")

# Step 3: Re-run the script from the virtual environment
if sys.executable != str(venv_python):
    print("🔁 Re-running inside virtual environment...")
    os.execv(str(venv_python), [str(venv_python), __file__])

# Step 4: Actual script starts here
print("🚀 Running script in virtual environment...")

# ... your original code continues here ...


import pandas as pd
from Bio.Seq import Seq
import time

def process_file(file_path, output_prefix):
    start_time = time.time()

    # Read the .tabular file
    df = pd.read_csv(file_path, sep='\t', header=None, names=['sequence'])

    # Extract UMI-Fw, UMI-Rv, and UMI pair
    df['UMI-Fw'] = df['sequence'].str[:15]
    df['UMI-Rv'] = df['sequence'].str[-15:]
    df['UMI pair'] = df['UMI-Fw'] + df['UMI-Rv']

    # Metrics before filtering
    reads_before = len(df)
    unique_umi_pair_before = df['UMI pair'].nunique()
    unique_umi_fw_before = df['UMI-Fw'].nunique()
    unique_umi_rv_before = df['UMI-Rv'].nunique()

    # Remove duplicates
    df.drop_duplicates(inplace=True)
    after_duplicates = len(df)

    # Extract DNA subsequence and translate
    df['extracted_sequence'] = df['sequence'].str[34:229]
    df['amino_acid_sequence'] = df['extracted_sequence'].apply(
        lambda x: str(Seq(x).translate()) if pd.notna(x) else ""
    )
    df['amino_acid_sequence'] = df['amino_acid_sequence'].str.replace(',', '', regex=False)

    # Extract CDR2 (first 10 aa) and CDR3 (last 16 aa)
    df['CDR2'] = df['amino_acid_sequence'].str[:10]
    df['CDR3'] = df['amino_acid_sequence'].str[-16:]

    # Remove stop codons
    df = df[~df['amino_acid_sequence'].str.contains(r'\*')]
    after_stop_codons = len(df)

    # Rearrange columns
    df = df[['UMI pair', 'UMI-Fw', 'UMI-Rv', 'sequence',
             'extracted_sequence', 'amino_acid_sequence',
             'CDR2', 'CDR3']]

    # Collect summary statistics for merged table
    summary_data = {
        "File": output_prefix,
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

    # Write counts for CDR2 and CDR3
    df['CDR2'].value_counts().reset_index(name='Count_CDR2').rename(columns={'index': 'CDR2'}) \
        .to_csv(f'{output_prefix}_cdr2_counts.csv', index=False)
    df['CDR3'].value_counts().reset_index(name='Count_CDR3').rename(columns={'index': 'CDR3'}) \
        .to_csv(f'{output_prefix}_cdr3_counts.csv', index=False)

    execution_time = time.time() - start_time
    print(f"Finished {file_path} in {execution_time:.2f} seconds")

    return summary_data


if __name__ == "__main__":
    summaries = []

    # Process n63.tabular to n85.tabular
    for i in range(63, 86):
        file_name = f"n{i}.tabular"
        if os.path.exists(file_name):
            summary = process_file(file_name, f"n{i}")
            summaries.append(summary)
        else:
            print(f"⚠️ Skipping {file_name} (not found)")

    # Merge summaries into one CSV
    if summaries:
        merged_df = pd.DataFrame(summaries)
        merged_df.to_csv("all_processing_summary.csv", index=False)
        print("✅ All summaries merged into all_processing_summary.csv")
    else:
        print("⚠️ No files processed, merged summary not created.")