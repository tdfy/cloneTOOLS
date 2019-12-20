library(tcR)
library(gridExtra)
library(grid)
library(circlize)
library(stringr)
library(dplyr)

path <- "c:/Export/TCRGD_2019Dec17/lib_1x2x3/"

#file_list <- c('1_S1.clonotypes.TRG.txt','2_S2.clonotypes.TRG.txt','3_S3.clonotypes.TRG.txt')

file_list <- list.files(path)

for (file in file_list){

  print(file)

  voodoo <- read.table(paste(path,file,sep=""), header=TRUE, sep="\t")  ##### Necessary to convert clone ID count to true integer***
  
  
  
  write.table(voodoo, file=paste(path,file,sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

}

#parse.file.list(file_list, 'mixcr')

emo <- parse.folder(path, 'mixcr')

name = "1x2x3"

png(paste(path,name,".png", sep=""))
prop <- vis.top.proportions(emo) 
penta_hope = paste(path,name, "_topCLONE.png", sep="")
ggsave(penta_hope)


twb.space <- clonal.space.homeostasis(emo)
clone <-vis.clonal.space(twb.space)
penta = paste(path,name, "_cloneSPACE.png", sep="")
ggsave(penta)

use <- vis.gene.usage(twb[[1]], HUMAN_TRGV, .main = 'Sample I V-usage')


stat <- cloneset.stats(emo)
write.table(stat, file=paste(path,name,"_stat.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

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

