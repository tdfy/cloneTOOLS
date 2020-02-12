library(tcR)
library(gridExtra)
library(grid)
library(circlize)
library(stringr)
library(dplyr)
library(VennDiagram)
library(rlist)
library(vwr)
lirbary(purrr)


path <- "c:/Export/TCRGD_2019Dec17/lib_289/convert"

setwd(path)

# file_list <- c('1_S1.clonotypes.TRG.txt','2_S2.clonotypes.TRG.txt','3_S3.clonotypes.TRG.txt')

file_list <- list.files(path)

for (file in file_list){

  voodoo <- read.table(paste(path,"/",file,sep=""), header=TRUE, sep="\t")  ##### Necessary to convert clone ID count to true integer***
  
  write.table(voodoo, file=paste(path,"/",file,sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

}

# emo <- parse.file.list(file_list, 'mixcr')

emo <- parse.folder(path, 'mixcr')

name = "2x8x9"

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
write.table(stat, file=paste(path,"/",name,"_stat.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

short_stat <- stat[,c("#Nucleotide clones","#Aminoacid clonotypes","#In-frames","#Out-of-frames","Mean.reads","Sum.reads")]
write.table(short_stat, file=paste(path,name,"_short_stat.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

##---------Top Clone Stats ------###

top1 <- emo[[1]][1:10,4:8]
write.table(top1, file=paste(path,name,"S2_gamma.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

top2 <- emo[[2]][1:10,4:8]
write.table(top2, file=paste(path,name,"S8_gamma.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

top3 <- emo[[3]][1:10,4:8]
write.table(top3, file=paste(path,name,"S9_gamma.txt",sep=""), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)

###------- Circos ----------------###


for (file in file_list){


  try <- file
  depth <- read.table(paste(path,"/",try,sep=""),header=T,sep="\t")
  
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
                           circos.text(mean(xlim), ylim[1], cex = 0.4, sector.name, facing = "clockwise", adj = c(-0.5, 1))
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
df1x <- as.data.frame(emo[1])
df1 <- as.data.frame(emo[1])
df1$Clone_ID <- rownames(df1)
colnames(df1)[6] <- "AA_Seq"

df2 <- as.data.frame(emo[2])
df2$Clone_ID <- rownames(df2)
colnames(df2)[6] <- "AA_Seq"

df3 <- as.data.frame(emo[3])
df3$Clone_ID <- rownames(df3)
colnames(df3)[6] <- "AA_Seq"


try3 <- sapply(df1$AA_Seq,as.factor)

try4 <- sapply(df2$AA_Seq,as.factor)

try5 <- sapply(df3$AA_Seq,as.factor)

ONE <- levels(try3)
TWO <- levels(try4)
THREE <- levels(try5)

namelist <- c("S2", "S8", "S9")

mix <-venn.plot <- draw.triple.venn(
  area1 = length(ONE),
  area2 = length(TWO),
  area3 = length(THREE),
  n12 = length(intersect(ONE,TWO)),
  n23 = length(intersect(TWO,THREE)),
  n13 = length(intersect(ONE,THREE)),
  n123 = length(Reduce(intersect, list(ONE,TWO,THREE))),
  category = namelist,
  fill = c("blue", "red", "green")) #,

gg <- grid.arrange(gTree(children=mix), top="TCR Gamma AA Clonotypes")

inter1 <- Reduce(intersect, list(ONE,TWO))
inter2 <- Reduce(intersect, list(ONE,THREE))
inter3 <- Reduce(intersect, list(TWO,THREE))
interALL <- Reduce(intersect,list(ONE,TWO,THREE))


top_share <- function (x,y) {
  # int_df_list<<- list()
  cols <- c(4:8,17)
  lib_name <- deparse(substitute(x))
  LN <<- paste("lib",lib_name,sep="_")

  nn1 <-str_replace_all(LN,"[\\[]","")
  LN <- str_replace_all(nn1,"[\\]]","")
  
  intersection <- deparse(substitute(y))
  interx1 <<- paste(intersection)
  n1 <- str_replace_all(interx1,"[\\[]","")
  interx <-str_replace_all(n1,"[\\]]","")
  

  sub_df <<- subset(x, x$AA_Seq %in% y)
  lib_flat <- apply(sub_df,2,as.character)
  
  int_name <<-paste(LN,interx,"all_share.txt",sep=".")
  # hope <<-list.append(int_df_list,int_name)
  
  write.table(lib_flat, file=paste(LN,interx,"all_share.txt",sep="."), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)
  
  
  print(paste(LN,interx,"sum",sep="."))
  
  newlib <<- aggregate(sub_df[,c(4)] ~ sub_df[,c(6)], data=sub_df, FUN=sum)
  assign(paste(LN,interx,"sum",sep="."), round(sum(newlib[2]), digits=3),envir = .GlobalEnv)

  # lib_flat <- apply(sub_df,2,as.character)
  lib_flatter <- lib_flat[1:10,cols]

  write.table(lib_flatter, file=paste(LN,interx,"top_share.txt",sep="."), sep="\t",row.names=FALSE,quote = FALSE,col.names=TRUE)
  
  png(paste(LN,interx,"top_share.png",sep="."), height=1000, width=2000)
  p<-tableGrob(lib_flatter)
  grid.arrange(p)
  dev.off()
}

df_list <- list(df1,df2,df1,df3,df2,df3)
intersX <- list(inter1,inter1,inter2,inter2,inter3,inter3)

mapply(top_share,df_list,intersX)




neo = matrix(c(100,lib_dots1L1L.dots2L1L.sum,lib_dots1L3L.dots2L3L.sum,lib_dots1L2L.dots2L2L.sum,100,lib_dots1L5L.dots2L5L.sum,lib_dots1L4L.dots2L4L.sum,lib_dots1L6L.dots2L6L.sum,100),nrow=3,ncol=3)
colnames(neo) <-paste(namelist,sep=",")
rownames(neo) <-paste(namelist,sep=",")

trip <- apply(neo,2,as.character)
write.table(trip, file='prop_matrix.txt', sep="\t",row.names=TRUE,quote = FALSE,col.names=TRUE)


##_____all intersections_______##

lib1 <- subset(df1, df1$X2_S2.CDR3.amino.acid.sequence %in% interALL)
newlib1 <- aggregate(lib1[,c(4)] ~ lib1[,c(6)], data=lib1, FUN=sum)
sum1all <- round(sum(newlib1[2]), digits=3)

lib2 <- subset(df2, df2$X8_S8.CDR3.amino.acid.sequence %in% interALL)
newlib2 <- aggregate(lib2[,c(4)] ~ lib2[,c(6)], data=lib2, FUN=sum)
sum2all <- round(sum(newlib2[2]), digits=3)


lib3 <- subset(df3, df3$X9_S9.CDR3.amino.acid.sequence %in% interALL)
newlib3 <- aggregate(lib3[,c(4)] ~ lib3[,c(6)], data=lib3, FUN=sum)
sum3all <- round(sum(newlib3[2]), digits=3)

##_______clone finder_____________##

df = NULL

leven_find <- function(x,y){
  # df = NULL
  
  xc <- cross2(x,twb[[y]]$CDR3.amino.acid.sequence)

  for (i in 1:length(xc)){
    score <- levenshtein.distance(xc[[i]][[1]],xc[[i]][[2]])
    if (score <5){
    Lev_Score <- c(score[[1]])
    MR1_AA <- c(xc[[i]][[1]])
    Clonotype <- c(xc[[i]][[2]])
    MR1_len <- nchar(xc[[i]][[1]])
    Clone_len <- nchar(xc[[i]][[2]])
    Adj_Lev <-score - abs(nchar(xc[[i]][[1]])-nchar(xc[[i]][[2]]))
    x <- data.frame(MR1_AA,Lev_Score,Clonotype,MR1_len,Clone_len, Adj_Lev)
    df<<-rbind(df,x)}else{}
  }
  return(df)
}

write.table(df, file=paste(path,"/","S6_delta_V35_MR1.csv",sep=""), sep=",",row.names=FALSE,quote = FALSE,col.names=TRUE)


