#### Set directory of script to working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Libraries:
for(i in 1){
  library(devtools)
  library(ggplot2)
  library(tidyr)
  library(grid)
  library(ggseqlogo)
  library(gridExtra)
  library(plyr)
  library(dplyr)
  library(magrittr)
  library(lattice)
  library("zoo")
  library(reshape2)
  library(stringr)
  library(ggsignif)
  library(ggpubr)
  library(xtable)
  library(corrplot)
  library(ggrepel)
  library(reticulate)
  library(rjson)
  library(RColorBrewer)
  library(seqinr)
  #install_github("vqv/ggbiplot")
  library(ggbiplot)
  library(readr)
  library(Biostrings)
  library(useful)
  library(Peptides)
  library(lemon)
  library(data.table)
}

#### Create Functions
for(i in 1){
  proportion <- function(x){
    rs <- sum(x);
    return(x / rs);
  }
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  substrVectorString <- function(x, y){
    characters<- c()
    for(i in 1:length(y)){
      aa <- substr(x, y[i],y[i])
      characters[i] <- aa
    }
    stringTotal <- paste(characters, collapse = "")
    return(stringTotal)
  }
  substrVectorColumn <- function(x, y){
    working <- as.data.frame(x)
    colnames(working) <- c("Seq")
    for(i in 1:length(y)){
      colname <- as.character(y[i])
      working[,colname] <- substr(working$Seq, y[i],y[i])
    }
    working$stringTotal <- do.call(paste0, working[as.character(y)])
    return(working$stringTotal)
  }
  split_str_by_index <- function(target, index) {
    index <- sort(index)
    substr(rep(target, length(index) + 1),
           start = c(1, index),
           stop = c(index -1, nchar(target)))
  }
  interleave <- function(v1,v2){
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
  }
  insert_str <- function(target, insert, index) {
    insert <- insert[order(index)]
    index <- sort(index)
    paste(interleave(split_str_by_index(target, index), insert), collapse="")
  }
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) {
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, nrow = ceiling(numPlots/cols))}
    if (numPlots==1) {print(plots[[1]])} else { grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      for (i in 1:numPlots) {
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }}}
  shift <- function(x, n){
    c(x[-(seq(n))], rep(0, n))
  }
  fudgeit <- function(){
    xm <- get("xm", envir = parent.frame(1))
    ym <- get("ym", envir = parent.frame(1))
    z  <- get("dens", envir = parent.frame(1))
    colramp <- get("colramp", parent.frame(1))
    image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =T)
  }
  Lab.palette <- colorRampPalette(c("white", "blue","yellow","red"), space = "Lab")
  MakeDensPlotsByFactor <- function(dataset, xvarstr, yvarstr, splitfactor, colspdf, xlimsvect, ylimsvect, densmax, xlabstr, ylabstr, mainstr, pdfname){
    levelsvector <- c(levels(dataset[,splitfactor]))
    print(levelsvector)
    nplots <- length(levelsvector)
    nrowspdf <- ceiling(nplots/colspdf)
    pdf(file=pdfname,height=(6*nrowspdf),width=(8*colspdf),onefile=TRUE,bg="white",version="1.7",useDingbats=FALSE)
    plots <- list()
    for(i in 1:nplots){
      sub <- dataset[dataset[,splitfactor]==levelsvector[i],]
      p1<-ggplot(data=sub,aes_string(xvarstr,yvarstr)) + 
        stat_density2d(aes(fill = ..density..^0.25), geom="tile", contour=FALSE) + 
        scale_fill_gradientn(colours=Lab.palette(20),guide="colorbar",limits=c(0,densmax))+
        geom_vline(xintercept = 0)+
        theme_bw() +
        theme(panel.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank() ,panel.border = element_rect(colour= "black"))+
        theme(legend.title = element_blank()) +
        theme(legend.key.height=(unit(2, "cm"))) +
        scale_x_continuous(limits=xlimsvect, expand = c(0, 0))+
        scale_y_continuous(limits=ylimsvect, expand = c(0, 0))+
        xlab(xlabstr) +
        ylab(ylabstr) +
        ggtitle(levelsvector[i])
      plots[[i]] <- p1
    }
    multiplot(plotlist=plots,cols=colspdf)
    dev.off()
  }
  MakeDensPlot <- function(dataset, xvarstr, yvarstr, xlimsvect, ylimsvect, densmax, xlabstr, ylabstr, mainstr,pdfname){
    p1<-ggplot(data=dataset,aes_string(xvarstr,yvarstr)) + 
      stat_density2d(aes(fill = ..density..^0.25), geom="tile", contour=FALSE) + 
      scale_fill_gradientn(colours=Lab.palette(20),guide="colorbar",limits=c(0,densmax))+
      geom_vline(xintercept = 0)+
      theme_bw() +
      theme(panel.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank() ,panel.border = element_rect(colour= "black"))+
      theme(legend.title = element_blank()) +
      theme(legend.key.height=(unit(2, "cm")))+
      scale_x_continuous(limits=xlimsvect, expand = c(0, 0))+
      scale_y_continuous(limits=ylimsvect, expand = c(0, 0))+
      xlab(xlabstr) +
      ylab(ylabstr) +
      ggtitle(mainstr)
    print(p1)
    ggsave(filename=pdfname, plot=p1, width=8, height=6)
    dev.off()
  }
  extract.with.context <- function(x, rows, after = 0, before = 0) {
    
    match.idx  <- which(rownames(x) %in% rows)
    span       <- seq(from = -before, to = after)
    extend.idx <- c(outer(match.idx, span, `+`))
    extend.idx <- Filter(function(i) i > 0 & i <= nrow(x), extend.idx)
    extend.idx <- sort(extend.idx)
    
    return(x[extend.idx, , drop = FALSE])
  }
  "%!in%" <- function(x,y)!("%in%"(x,y))
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  # Median test
  median.test <- function(x, y){
    z <- c(x, y)
    g <- rep(1:2, c(length(x), length(y)))
    m <- median(z)
    fisher.test(z < m, g)$p.value
  }
  zeroVar <- function(data, useNA = "ifany") {
    out <- apply(data, 2, function(x) {length(table(x, useNA = useNA))})
    which(out==1)
  }
  scaleFUN <- function(x) sprintf("%.1f", x)
  aaCheck <- function(seq){
    if(!any(lengths(seq) > 1)){
      seq <- toupper(seq)
      seq <- gsub(pattern = "[[:space:]]+",replacement = "",x = seq)
      seq <- strsplit(x = seq,split = "")
    } else {
      seq <- lapply(seq,toupper)
    }
    check <- unlist(lapply(seq,function(sequence){
      !all(sequence%in%c("A" ,"C" ,"D" ,"E" ,"F" ,"G" ,"H" ,"I" ,"K" ,"L" ,"M" ,"N" ,"P" ,"Q" ,"R" ,"S" ,"T" ,"V" ,"W" ,"Y", "-"))
    }))
    if(sum(check) > 0){
      sapply(which(check == TRUE),function(sequence){warning(paste0("Sequence ",sequence," has unrecognized amino acid types. Output value might be wrong calculated"),call. = FALSE)})
    }
    return(seq)
  }
  aaComp2 <- function (seq) 
  {
    seq <- aaCheck(seq)
    seq <- lapply(seq, function(seq) {
      table(unlist(seq))
    })
    aacomp <- lapply(seq, function(seq) {
      AA <- matrix(ncol = 2, nrow = 9)
      rownames(AA) <- c("Tiny", "Small", "Aliphatic", "Aromatic", 
                        "NonPolar", "Polar", "Charged", "Basic", "Acidic")
      colnames(AA) <- c(".n", ".p")
      AA[1, 1] <- sum(seq[c("A", "C", "G", "S", "T")], na.rm = TRUE)
      AA[2, 1] <- sum(seq[c("A", "B", "C", "D", "G", "N", "P", 
                            "S", "T", "V")], na.rm = TRUE)
      AA[3, 1] <- sum(seq[c("A", "I", "L", "V")], na.rm = TRUE)
      AA[4, 1] <- sum(seq[c("F", "H", "W", "Y")], na.rm = TRUE)
      AA[5, 1] <- sum(seq[c("A", "C", "F", "G", "I", "L", "M", 
                            "P", "V", "W", "Y")], na.rm = TRUE)
      AA[6, 1] <- sum(seq[c("D", "E", "H", "K", "N", "Q", "R", 
                            "S", "T", "Z")], na.rm = TRUE)
      AA[7, 1] <- sum(seq[c("B", "D", "E", "H", "K", "R", "Z")], 
                      na.rm = TRUE)
      AA[8, 1] <- sum(seq[c("H", "K", "R")], na.rm = TRUE)
      AA[9, 1] <- sum(seq[c("B", "D", "E", "Z")], na.rm = TRUE)
      AA[, 2] <- (AA[, 1]/sum(seq) * 100)
      AA <- round(AA, 3)
      AA <- melt(t(AA))
      AA$Name1 <- paste(AA$Var2,AA$Var1,sep="")
      AA <- AA[,c("Name1","value")]
      AA <- setNames(data.frame(t(AA[,-1])), AA[,1])
      return(AA)
    })
    
    return(aacomp)
  }
}

cols <- c("#F8766D", "#619CFF", "#00BA38")

Vsig4_CDR23_LibData <- read.delim("Vsig4_CDR23_Library_data.txt",stringsAsFactors = T,header=T)
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

Vsig4_CDR23_LibData_Filtered <- subset(Vsig4_CDR23_LibData, dG_ac < (-17.5) & IntraClash_binder < (10) &  TANGO_Mut_N < 300 & CavityVolumeAb < 75 & CavityVolumeInt < 125)

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

