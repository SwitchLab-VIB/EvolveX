#### Set directory of script to working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Libraries:
for(i in 1){
  library(ggplot2)
  library(tidyr)
  library(openxlsx)
  library(plyr)
  library(dplyr)
  library(ggpubr)
}

cols <- c("#F8766D", "#619CFF", "#00BA38")

Vsig4_CDR23_LibData <- read.xlsx("InSilicoVariables.xlsx",sheet="Sheet1")
Vsig4_CDR23_LibData <- subset(Vsig4_CDR23_LibData,goal=="AffinityMouse" & Type == "Candidates_Trp101Fix" & Backbone == "RP5imoBA")
HitSeqs <- c("QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYGMGWFRQAPGKEREFVAAIRWNGVETYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGRWDLFGSYDQDNYEYWGQGTQVTVSS",
             "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYGMGWFRQAPGKEREFVAAIRWNGVETYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGRWDLFGSFDQDNYEYWGQGTQVTVSS")
Vsig4_CDR23_LibData$HitMouse <- "N"
Vsig4_CDR23_LibData$HitMouse[Vsig4_CDR23_LibData$MutSeq_N %in% HitSeqs] <- "Y"
Vsig4_CDR23_LibData <- Vsig4_CDR23_LibData %>% mutate(dG_rank = order(order(dG, decreasing=FALSE)))
Vsig4_CDR23_LibData <- Vsig4_CDR23_LibData %>% mutate(dGac_rank = order(order(dG_ac, decreasing=FALSE)))
Vsig4_CDR23_LibData <- Vsig4_CDR23_LibData %>% mutate(dGbinder_rank = order(order(dG_binder, decreasing=FALSE)))
Vsig4_CDR23_LibData <- Vsig4_CDR23_LibData %>% mutate(IntraClashes_rank = order(order(IntraClash_binder, decreasing=FALSE)))
Vsig4_CDR23_LibData <- Vsig4_CDR23_LibData %>% mutate(TANGO_rank = order(order(TANGO_Mut_N, decreasing=FALSE)))
Vsig4_CDR23_LibData$RankSum <- Vsig4_CDR23_LibData$dGac_rank+Vsig4_CDR23_LibData$dGbinder_rank+Vsig4_CDR23_LibData$IntraClashes_rank
Vsig4_CDR23_LibData <- Vsig4_CDR23_LibData %>% mutate(RankTotal = order(order(RankSum, decreasing=FALSE)))

for(i in 1){
  HitSeqs_Offrate_LibData <- subset(Vsig4_CDR23_LibData,MutSeq_N %in% HitSeqs)
  HitSeqs_Offrate_LibData$Name <- mapvalues(HitSeqs_Offrate_LibData$Name,from = c("Rob_6731","Rob_6732"),to=c("VHH_m4","VHH_m2"))
  HitSeqs_Offrate_LibData$Name <- factor(HitSeqs_Offrate_LibData$Name,levels=c("VHH_m2","VHH_m4"))
  Library_LibData <- Vsig4_CDR23_LibData
  colnames_toPlot <- colnames(select_if(Vsig4_CDR23_LibData, is.numeric))
  colnames_toExclude <- c("ThreadName","NextIteration","CavityVolumeTotal","CavityVolumeLig","ChargeLig","PiPiNum","PiPiDist",
                          "IntraclashesGroup2","sloop_entropy","mloop_entropy","cis_bond","helix.dipole","water.bridge","disulfide","partial.covalent.bonds","Entropy.Complex","Number.of.Residues",
                          "Interface.Residues.BB.Clashing","ThreadName2","Dum2","CavityVolumeInt",
                          "ddG","ddG_binder","dIntraClash_binder","ddG_ac","TANGO_WT_N","dTANGO_Mut_N","SeqLen","MutNameLen","DecisionNum","SeqSim80","SeqSimNorm80","SeqID","SeqSimSelf80",
                          "Interface.Residues.Clashing","Interface.Residues.VdW.Clashing","energy.Ionisation","electrostatic.kon","backbone.clash","torsional.clash",
                          "Iteration","PiPiEn","Ionic","Hydrophobic","ZS_CavAb","ZS_CavInt","ZS_ChAb","ZS_Sum","IntraclashesGroup1","Solvation.Polar","Solvation.Hydrophobic","entropy.sidechain","entropy.mainchain","Interface.Residues","IntraClash_binder",
                          "dG_rank","dGac_rank","dGbinder_rank","IntraClashes_rank","TANGO_rank","RankSum","RankTotal","InterfaceSurfaceArea","dG_ac")
  colnames_toPlot <- colnames_toPlot[! colnames_toPlot %in% colnames_toExclude]
  colnames_toPlot <- c("Interaction.Energy","CationPi","Backbone.Hbond","Sidechain.Hbond","Electrostatics","Van.der.Waals","Van.der.Waals.clashes","dG","dG_binder","TANGO_Mut_N","CavityVolumeAb","ChargeAb")
  colnames_labels <- c("Interaction DG (kcal/mol)","Interaction Cation-Pi","Interaction Backbone H-bond (kcal/mol)","Interaction Sidechain H-bond (kcal/mol)","Interaction Electrostatics (kcal/mol)","Interaction VanDerWaals (kcal/mol)","Interaction VanDerWaals Clashes (kcal/mol)","Complex DG (kcal/mol)","VHH DG (kcal/mol)","VHH TANGO","VVH Cavity Volume","VHH Net Charge")
  plotlist_LibData <- list()
  for(j in 1:length(colnames_toPlot)){
    df_all_temp <- Library_LibData[,c("Name",colnames_toPlot[j])]
    df_all_toPlot <- df_all_temp
    colnames(df_all_toPlot) <- c("Name","Variable")
    df_hits_temp <- HitSeqs_Offrate_LibData[,c("Name",colnames_toPlot[j])]
    df_hits_toPlot <- df_hits_temp
    colnames(df_hits_toPlot) <- c("Name","Variable")
    p1 <- ggplot(df_all_toPlot, aes(x=Variable)) +
      geom_density()+
      theme(legend.background = element_blank(),legend.position = "top",legend.box.background = element_blank(),legend.title=element_blank(),legend.key = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour="black"))+
      # geom_vline(data=df_hits_toPlot, aes(xintercept=Variable, color=Name),
      #            linetype="dashed")+
      geom_segment(data = df_hits_toPlot,
        mapping = aes(x = Variable, xend = Variable,y = -Inf, yend = Inf,
          color = factor(Name),linetype = factor(Name)))+
      scale_linetype_manual(values = c("11", "33")) +  # two different dash patterns
      labs(title="",x=colnames_labels[j],y="Density")
    plotlist_LibData[[j]] <- p1
  }
  pdf(file=paste("HitLocation.pdf",sep=""),height=10,width=12,onefile=TRUE,bg="white",version="1.7",useDingbats=FALSE)
  plots <- ggarrange(plotlist=plotlist_LibData,
                     labels = LETTERS[1:length(plotlist_LibData)],
                     ncol = 3, nrow = 4,common.legend = TRUE, legend="top")
  print(plots)
  dev.off()
}

