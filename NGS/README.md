# Vsig4-CDR23-PCR-based-NGS-FoldX-analysis
NGS and FoldX analysis of phage selected Vsig4 targeting in silico design VHH sequences.

VHH CDR2 and CDR3 sequences were designed based on existing VHH_WT bound to both mouse and human Vsig4 using EvolveX pipeline. The phage library with designed sequences was selected with several rounds of mouse or human Vsig4. Traditional PCR-prep Illumina NGS was performed for the input and output from the selection campaigns. Sequences were further assessed with FoldX. Here are the analysis protocol and scripts we used.

NGS raw data were preliminarily processed at Galaxy website with PEAR merger. The sequence column was cut and the sequence only data were save as nxx.tabular (for example: n66.tabular). Sample information was stored in **Sample_infor.xlsx**.

Then these files were further processed with in-house python script **count.py**. UMI information was used to deduplicate and CDR2, CDR3 sequences were extracted. For each nxx.tabular file, there were two output csv files: one contained unique CDR2 sequences and their counts; one contained unique CDR3 sequences and their counts. There was also a summarized processing file for all the samples. Output files were stored in the folder **processed_file_after_count.zip**. 

Afterward, files with CDR sequence and counts can be further merged with their negative control groups, sorted and mapped based on customized requirement (**file_list.csv**) and designed sequence map (**id.csv**).Python script **merge.py** was used. Output files were stored in the folder **processed_file_after_merge.zip**. 

Plots for count and enrichment for each CDR sequence in each round were made in the folder **Plots**. 

Figures for sequence logos were made in the folder **Logo**. 

Further FoldX modeling and anaylysis were made in the folder **FoldX_analysis**.
