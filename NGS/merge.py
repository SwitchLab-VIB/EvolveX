import csv

def merge_filter_and_map(cdr_counts_file1, cdr_counts_file2, id_file, output_file,
                         seq_type="CDR2", min_count=1, min_enrichment=1):
    """
    Merge counts from two files, filter (> min_count AND > min_enrichment), sort, and map IDs from id.csv.
    
    seq_type: "CDR2" or "CDR3"
    """
    # --- Step 1: read both count files ---
    data1, data2 = {}, {}
    for file, data in [(cdr_counts_file1, data1), (cdr_counts_file2, data2)]:
        with open(file, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # skip header
            for row in reader:
                sequence = row[0]
                count = int(row[1])
                data[sequence] = count

    # --- Step 2: merge + calculate enrichment ---
    merged_data = {}
    for seq in set(data1.keys()) | set(data2.keys()):
        count1 = data1.get(seq, 0)
        count2 = data2.get(seq, 0)
        enrichment = count1 / count2 if count2 != 0 else count1 / 1
        merged_data[seq] = [count1, count2, enrichment]

    # --- Step 3: filter + sort ---
    rows = []
    for seq, (count1, count2, enrichment) in merged_data.items():
        if count1 > min_count and enrichment > min_enrichment:  # ✅ AND filter
            rows.append({
                'Sequence': seq,
                'Count_binding': count1,
                'Count_nsbinding': count2,
                'Enrichment_factor': enrichment
            })
    rows.sort(key=lambda x: x['Count_binding'], reverse=True)

    # --- Step 4: load ID mappings (keep first match only) ---
    id_map = {}
    with open(id_file, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            cdr2, cdr3, identifier = row
            if seq_type == "CDR2":
                if cdr2 not in id_map:   # ✅ keep only the first match
                    id_map[cdr2] = identifier
            elif seq_type == "CDR3":
                if cdr3 not in id_map:   # ✅ keep only the first match
                    id_map[cdr3] = identifier

    # --- Step 5: map IDs ---
    for row in rows:
        seq = row['Sequence']
        row['ID'] = id_map.get(seq, "")

    # --- Step 6: write final file ---
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Sequence', 'Count_binding', 'Count_nsbinding', 
                                               'Enrichment_factor', 'ID'])
        writer.writeheader()
        writer.writerows(rows)

    print(f"✅ Final file saved: {output_file}")


def run_batch(file_list_csv, id_file):
    """
    Read the file list CSV and process each row automatically for CDR2 and CDR3 files.
    file_list_csv: CSV with columns -> file1,file2,output,min_count,min_enrichment
    """
    with open(file_list_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            file1_cdr2 = f"{row['file1']}_cdr2_counts.csv"
            file2_cdr2 = f"{row['file2']}_cdr2_counts.csv"
            output_cdr2 = f"{row['output']}_CDR2.csv"
            min_count = int(row['min_count'])
            min_enrichment = int(row['min_enrichment'])
            merge_filter_and_map(file1_cdr2, file2_cdr2, id_file, output_cdr2,
                                 seq_type="CDR2", min_count=min_count, min_enrichment=min_enrichment)

            file1_cdr3 = f"{row['file1']}_cdr3_counts.csv"
            file2_cdr3 = f"{row['file2']}_cdr3_counts.csv"
            output_cdr3 = f"{row['output']}_CDR3.csv"
            merge_filter_and_map(file1_cdr3, file2_cdr3, id_file, output_cdr3,
                                 seq_type="CDR3", min_count=min_count, min_enrichment=min_enrichment)


# --- Example usage ---
if __name__ == "__main__":
    run_batch("file_list.csv", "id.csv")
