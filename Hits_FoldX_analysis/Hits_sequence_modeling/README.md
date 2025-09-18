# Hits re-modeling from new Ala-starting points

Put your Ala_starting-point pdb structures generated in the folder **Hits_FoldX_analysis/Ala_starting** (5imkBA_Repair_1.pdb, 5imlBA_Repair_1.pdb, 5immBA_Repair_1.pdb, 5imoBA_Repair_1.pdb, 9ezh_Repair_1.pdb, 9ezi_Repair_1.pdb, 9ezu_Repair_1.pdb, 9ezv_Repair_1.pdb, 9ezw_Repair_1.pdb) in this folder, then run the bash file **Hits_BM.sh** first to generate the complex structure models for all the hits based on different Ala-starting points (backbones). Three models will be generated for each hit each Ala-starting point, so the average energies can be calculated manually afterward. 

Then run bash file **Hits_AC.sh** to calculate the energies of all generated models. Summary files (e.g. Summary_5imk_AC.fxout) contained the VHH individual energies and interaction energies. The average energies were then calculated manually. 

