library(tcR)
library(gridExtra)
library(grid)
library(circlize)
library(stringr)
library(dplyr)
library(VennDiagram)


path <- "c:/Export/TCRGD_2019Dec17/lib_1x2x3"

setwd(path)

file_list <- c('1_S1.clonotypes.TRG.txt','2_S2.clonotypes.TRG.txt','3_S3.clonotypes.TRG.txt')

file_list <- list.files(path)

for (file in file_list){

  print(file)

  voodoo <- read.table(paste(path,file,sep=""), header=TRUE, sep="\t")  ##### Necessary to convert clone ID count to true integer***
  
  
  
  write.table(voodoo, file=paste(path,file,sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

}

emo <- parse.file.list(file_list, 'mixcr')

emo <- parse.folder(path, 'mixcr')

name = "1x2x3"

prop <- vis.top.proportions(emo) 
penta_hope = paste(path,name, "_topCLONE.png", sep="")
ggsave(penta_hope)


twb.space <- clonal.space.homeostasis(emo)
clone <-vis.clonal.space(twb.space)
penta = paste(path,name, "_cloneSPACE.png", sep="")
ggsave(penta)

png("lib_gen_v.png", width = 1000, height = 1000)
v <- vis.gene.usage(emo, HUMAN_TRGV, .main = 'Variable Gene Distribution', .dodge = F, .ncol = 3)
plot(v)
dev.off()
grid.newpage();

png("lib_gen_j.png", width = 1000, height = 1000)
j <-vis.gene.usage(emo, HUMAN_TRGJ, .main = 'Joining Gene Distribution', .dodge = F, .ncol = 3)
plot(j)
dev.off()
grid.newpage();


stat <- cloneset.stats(emo)
write.table(stat, file=paste(path,name,"_stat.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

short_stat <- stat[,c("#Nucleotide clones","#Aminoacid clonotypes","#In-frames","#Out-of-frames","Mean.reads","Sum.reads")]
write.table(short_stat, file=paste(path,name,"_short_stat.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

##---------Top Clone Stats ------###

top1 <- emo[[1]][1:10,4:8]
write.table(top1, file=paste(path,name,"S1_gamma.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

top2 <- emo[[2]][1:10,4:8]
write.table(top2, file=paste(path,name,"S2_gamma.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

top3 <- emo[[3]][1:10,4:8]
write.table(top3, file=paste(path,name,"S3_gamma.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

###------- Circos ----------------###


for (file in file_list){


  try <- file
  depth <- read.table(paste(path,try,sep=""),header=T,sep="\t")
  
  depth$vgene <- str_split_fixed(depth$allVHitsWithScore,"[*]",4)[,1]
  depth$jgene <- str_split_fixed(depth$allJHitsWithScore,"[*]",4)[,1]
  
  newdf <- depth
  newdf <-select(newdf,cloneFraction,vgene,jgene)
  
  agg <- aggregate(cloneFraction ~ vgene+jgene, newdf, FUN=sum)
  
  chordDiagram(agg,annotationTrack = c("axis","grid"), preAllocateTracks = list(track.height = 0.2), transparency = 0.5)
  
  circos.trackPlotRegion(track.index = 1, bg.border = NA,
                         panel.fun = function(x, y) {
                           sector.name = get.cell.meta.data("sector.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(mean(xlim), ylim[1], cex = 0.65, sector.name, facing = "clockwise", adj = c(-0.5, 1))
                         }
  )
  
  for(si in get.all.sector.index()) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
                facing = "bending.inside", niceFacing = TRUE, col = "white",adj=c(3,-5))
  }
  
  chord_plot <-title(paste(file))
  
  dev.copy2pdf(file=paste(path,file,".pdf",sep=""))

}


###----------------- Venn ---------------------###
df1 <- as.data.frame(emo[1])
df2 <- as.data.frame(emo[2])
df3 <- as.data.frame(emo[3])

try3 <- sapply(df1$X1_S1.CDR3.amino.acid.sequence,as.factor)

try4 <- sapply(df2$X2_S2.CDR3.amino.acid.sequence,as.factor)

try5 <- sapply(df3$X3_S3.CDR3.amino.acid.sequence,as.factor)

ONE <- levels(try3)
TWO <- levels(try4)
THREE <- levels(try5)

mix <-venn.plot <- draw.triple.venn(
  area1 = length(ONE),
  area2 = length(TWO),
  area3 = length(THREE),
  n12 = length(intersect(ONE,TWO)),
  n23 = length(intersect(TWO,THREE)),
  n13 = length(intersect(ONE,THREE)),
  n123 = length(Reduce(intersect, list(ONE,TWO,THREE))),
  category = c("S1", "S2", "S3"),
  fill = c("blue", "red", "green")) #,

gg <- grid.arrange(gTree(children=mix), top="TCR Gamma AA Clonotypes")

inter <- Reduce(intersect, list(ONE,THREE))

threeish <- subset(df1, df1$X1_S1.CDR3.amino.acid.sequence %in% inter)
sum(threeish$X1_S1.Read.proportion)

fourish <- subset(df2, df2$X3_S2.CDR2.amino.acid.sequence %in% inter)
sum(fourish$X2_S2.Read.proportion)

fiveish <- subset(df3, df3$X3_S3.CDR3.amino.acid.sequence %in% inter)
sum(fiveish$X3_S3.Read.proportion)

new3 <- select(threeish,X1_S1.Read.proportion,X1_S1.CDR3.amino.acid.sequence)
new3a <- aggregate(X1_S1.Read.proportion ~ X1_S1.CDR3.amino.acid.sequence, data=new3, FUN=sum)
