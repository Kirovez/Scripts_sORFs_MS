
#setwd("C:/Users/Georgij/Downloads/")

setwd("D:/Rworking/MSspectra")
#.libPaths("D:\\Rpackages")
#########   Upload data about identified peptide sequences   #########
#install.packages("xlsx")
#########   Upload data from experiment mgf-file   #########
#source("https://bioconductor.org/biocLite.R")
#biocLite("MSnbase")
library(xlsx)
library(ggplot2)
library(MSnbase)
library(ggthemes)
library("scales")
library(grid)
require("ggrepel")

singleMSpeptides = fread("D:\\PycharmProjects\\sORF\\singleMSpeptide.txt", header = F)

PeptIdent <- read.xlsx("D:\\Rworking\\MSspectra\\mgf data\\end_peptides_orbitrap.xlsx", 1, stringsAsFactors = FALSE)

all_pep = as.data.table(rbind(PeptIdent, ProttIdent))
count_pep = all_pep[,.(nMS = length(Raw.file)), by = "Sequence"]

count_pep[nMS != 1 & Sequence %in% singleMSpeptides$V1]

additional_single = count_pep[nMS == 1]

singleMSpeptides[!(V1 %in% additional_single$Sequence) & !(V1 %in% all_pep$Sequence)]

single_everywhere_orbitrap = intersect(additional_single$Sequence, singleMSpeptides$V1)

count_pep[Sequence %in% singleMSpeptides$V1]

count_pep = count_pep[Sequence %in% singleMSpeptides$V1]

nrow(PeptIdent)

PeptIdent = as.data.table(PeptIdent)[Sequence %in% count_pep$Sequence]

length(unique(PeptIdent$Sequence))

Mgf <- readMgfData("D:\\Rworking\\MSspectra\\mgf data\\end_peptides.txt_selectedMGF.mgf")
length(Mgf)

i = 5

#########   Plot spectra and export pdf-file   #########
pdf(file = "./SupplementaryFile4.pdf", width = 10, height = 7)
for(i in 1:nrow(PeptIdent)) {
  Seq <- PeptIdent$Sequence[i]
  Spec <- PeptIdent$Raw.file[i]
  scan <- PeptIdent$Scan.number[i]
  sorf =  PeptIdent$sORFs[i]
  pep = PeptIdent$PEP[i]
  score = PeptIdent$Score[i]
 
  
  if(length(which(grepl(Spec, fData(Mgf)[,1]) & grepl(scan, fData(Mgf)[,1]))) == 1 & Seq %in% single_everywhere_orbitrap) {
    PeakList <- Mgf[[which(grepl(Spec, fData(Mgf)[,1]) & grepl(scan, fData(Mgf)[,1]))]]
    
    
    # all peaks
    all=data.frame(mz = mz(PeakList), intensity = 
                     intensity(PeakList))
    
    # peaks used for identification
    cf = calculateFragments(Seq, PeakList, type=c("a", "b", "y"))
    
    ## choose colors
    cf$col = ifelse(grepl("y", cf$type), "blue","red")
    
    if (nrow(cf[grepl("a", cf$type),]) > 0){
      cf[grepl("a", cf$type),]$col = "green"
    }
    
    cf$col = as.factor(cf$col)
 
    # table with annotation for peaks
    ann = data.frame(x = cf$mz, y =cf$intensity, text = cf$ion) # type column can be used for labels
    ann$y = ann$y + mean(ann$y)*0.1 # add space between peak end and lablel. Should be improved as it is sometime overlapped with a peak!

    # plotting
    par(mar = c(5.1, 4.1, 1.1, 1.1))
    print(ggplot() + 
      geom_segment(aes(x=mz, xend=mz, y=0, yend=intensity), data=all, colour = "grey",size = 0.3) + # lines for ALL peaks. 
      geom_point(aes(x=mz, y=intensity, col = col), size = 1.5, data=cf)  + # point in the end of the peak. Note: it can be removed
      geom_segment(aes(x=mz, xend=mz, y=0, yend=intensity,col=col), data=cf, size = 0.4) + # lines for the only  peaks used for identification
      theme_classic() + 
        xlim(0,max(all$mz)) +
        ggtitle(
          paste(paste(paste("Peptide:", Seq), paste("sORF:", sorf), sep=" "), 
                paste(paste("PEP:", pep), paste("Score:", score), paste("rt:",rtime(PeakList)), paste("peaks:",peaksCount(PeakList)), sep = " "),
          sep = "\n")) +
      ylab("Intensity") + 
      xlab("m/z") + 
      theme(legend.position="none", text = element_text(size=20, family = "sans"),
            axis.text.x = element_text(size=15, hjust = 1, color = "black"),
            plot.title = element_text(colour = "darkblue", size=13)) + 
        geom_text_repel(data=cf, aes(mz,intensity,label=ion, col=col), vjust=1.6, size=3.5)  +# labels +
      scale_color_manual(values = levels(cf$col)) 
        
    )
  } else if (Seq %in% singleMSpeptides$V1) {
    print(Seq)
    if(length(which(gsub('"', "", fData(Mgf)[,1]) == Spec)) == 0){
      print(paste("No specrum for peptide:", Seq, sep = " "))
    
  } else if(length(which(gsub('"', "", fData(Mgf)[,1]) == Spec)) > 1){
    print(paste("Too many specra for peptide:", Seq, sep = " "))
  }
} 
    
}

dev.off()

remove(i, Seq, Spec, Mgf, PeakList, PeptIdent)

dat <- data.frame(x=runif(10),y=runif(10),
                  grp = rep(LETTERS[1:5],each = 2),stringsAsFactors = TRUE)
library(RColorBrewer)
myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(cf$col)
colScale <- scale_colour_manual(name = "grp",values = myColors)
