# NGS enriched / negative CDR2 sequence re-modeling from new Ala-starting points

Put your Ala_starting-point pdb structures generated in the folder **NGS/FoldX_analysis/Ala_starting/CDR2** (5imkBA_Repair_1.pdb, 5imlBA_Repair_1.pdb, 5immBA_Repair_1.pdb, 5imoBA_Repair_1.pdb, 9ezh_Repair_1.pdb, 9ezi_Repair_1.pdb, 9ezu_Repair_1.pdb, 9ezv_Repair_1.pdb, 9ezw_Repair_1.pdb) in this folder, then run the bash file **NGS_BM.sh** first to generate the complex structure models for all the NGS enriched / negative CDR2 sequences based on different Ala-starting points (backbones). 

For hVsig4 complexes (5imk, 5iml, 9ezh, 9ezi, 9ezu, 9ezv, 9ezw): #1-38 are models re-built based on NGS enriched CDR2 sequences, #39-76 are models re-built based on negative CDR2 sequences.
For mVsig4 complexes (5imm, 5imo): #1-31 are models re-built based on NGS enriched CDR2 sequences, #32-62 are models re-built based on negative CDR2 sequences.

Then run bash file **NGS_AC.sh** to calculate the energies of all generated models. Summary files (e.g. **Summary_5imk_cdr2_AC.fxout**) contained the VHH individual energies and interaction energies. The ddG for VHH individual energy and interaction energy were then calculated manually. 

All results ddG data was summarized in the excel file **NGS/FoldX_analysis/NGS_sequence_modeling/SUMMARY.xlsx**.
