# Hits_sequence/structure FoldX analysis

With available sequences and structures of hits in complex with hVsig4 or mVsig4, we did the following FoldX modeling and analysis:

Sidechains in all nine available structures (PDBID: 5IMM and 5IMO for VHH_mVsig4, and 5IMK, 5IML, 9EZU, 9EZV, 9EVW, 9EZH, and 9EZI for VHH_hVsig4) were initially optimized by using the FoldX command RepairPDB. Structure models of all VHH candidates including VHH_WT, in complex with their respective target were generated using the FoldX BuildModel function, by first mutating the explored residues to Alanine (Folder **Ala_starting**) and then mutating to the corresponding amino acids (Folder **Hits_sequence_modeling**). The interaction energy (ΔGinteraction) and VHH individual energy (ΔGVHH) were calculated using FoldX AnalyseComplex function

