# NGS enriched / negative CDR2 sequence re-modeling from new Ala-starting points

Put your Ala_starting-point pdb structures generated in the folder **NGS/FoldX_analysis/Ala_starting/CDR2** (5imkBA_Repair_1.pdb, 5imlBA_Repair_1.pdb, 5immBA_Repair_1.pdb, 5imoBA_Repair_1.pdb, 9ezh_Repair_1.pdb, 9ezi_Repair_1.pdb, 9ezu_Repair_1.pdb, 9ezv_Repair_1.pdb, 9ezw_Repair_1.pdb) in this folder, then run the bash file **NGS_BM.sh** first to generate the complex structure models for all the NGS enriched / negative CDR2 sequences based on different Ala-starting points (backbones). 

Then run bash file **NGS_AC.sh** to calculate the energies of all generated models. Summary files (e.g. **Summary_5imk_cdr2_AC.fxout**) contained the VHH individual energies and interaction energies. The ddG for VHH individual energy and interaction energy were then calculated manually. 

All results ddG data was summarized in the excel file **NGS/FoldX_analysis/NGS_sequence_modeling/SUMMARY.xlsx**.
