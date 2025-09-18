# NGS sequence modelling based on available structures

With available sequences and structures of hits in complex with hVsig4 or mVsig4, we did the following FoldX modeling and analysis:

Sidechains in all nine available structures (PDBID: 5IMM and 5IMO for VHH_mVsig4, and 5IMK, 5IML, 9EZU, 9EZV, 9EVW, 9EZH, and 9EZI for VHH_hVsig4) were initially optimized by using the FoldX command _RepairPDB_. Structure models incorporating enriched or negative CDR sequences were generated using the FoldX _BuildModel_ function, by first mutating the explored residues to Alanine (Folder **Ala_starting**) and then mutating to the corresponding amino acids (Folder **NGS_sequence_modeling**). The interaction energy change (ΔΔGinteraction) and VHH individual energy change (ΔΔGVHH) were calculated using FoldX _AnalyseComplex_ function. 

