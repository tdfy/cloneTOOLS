library(immunarch)
library(circlize)
library(stringr)
library(vegan)


file_path = "c:/Export/30-359090205/chain_loc/gamma/UPCC21413"

immdata <- repLoad(file_path)

exp_vol = repExplore(immdata$data, .method = "volume")
exp_cnt = repExplore(immdata$data, .method = "count")

imm_pr = repClonality(immdata$data, .method = "clonal.prop")

imm_top = repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top

exp_len = repExplore(immdata$data, .method = "len", .col = "aa") #Meta break down of CDR3 length distributions
p6 = vis(exp_len, .by=c("Method"), .meta=immdata$meta,.test = F)

div_insimp = repDiversity(immdata$data, "inv.simp") #diversity
p1 = vis(div_insimp)



# ##-----------------------------------------------------------------------------------------------##



imm_ov1 = repOverlap(immdata$data, .method = "public", .verbose = F, .col ="aa")
imm_ov2 = repOverlap(immdata$data, .method = "overlap", .verbose = F,.col ="aa")


imm_gu = geneUsage(immdata$data, "hs.trgv", .norm = T)

vis(imm_gu)


vis(imm_gu, .plot = "tree")

imm_gu_cor = geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F) #gene usage heat map



p1 = vis(spectratype(immdata$data[[1]], .quant = "count", .col = "aa+v"))

grid.arrange(p1, p2,p3,p4,p5,p6,p7,p8, ncol = 4)
#---#
imm_gu_J = geneUsage(immdata$data, "hs.trgj", .norm = T)
imm_gu_cor_J = geneUsageAnalysis(imm_gu_J, .method = "cor", .verbose = F) #gene usage FOR JOING SEGMENT
vis(imm_gu_cor_J)


# ##-----------------------------------------------------------------------------------------------##


pr.aav = pubRep(immdata$data[c(1,5)], "aa+v", .verbose = F)



# ##-----------------------------------------------------------------------------------------------##

df1 <- as.data.frame(immdata$data[[3]])
df2 <- as.data.frame(immdata$data[[4]])

df3 <- as.data.frame(immdata$data[[5]])
df3 <- as.data.frame(immdata$data[[5]])



try3 <- levels(sapply(immdata$data[[3]]$CDR3.aa,as.factor))

try4 <- levels(sapply(immdata$data[[4]]$CDR3.aa,as.factor))

try5 <- levels(sapply(immdata$data[[5]]$CDR3.aa,as.factor))

try6 <- levels(sapply(immdata$data[[1]]$CDR3.aa,as.factor))



public <- intersect(intersect(try3,try4),intersect(try5,try6))


three <- list(immdata$data[[3]]$CDR3.aa)
four <- list(immdata$data[[4]]$CDR3.aa)
five <- list(immdata$data[[5]]$CDR3.aa)



##-----------------------------------------------------------------------##

target = immdata$data[[5]] %>% select(CDR3.aa, V.name) %>% head(20) #top clone target option

target = public[1:20]
tc = trackClonotypes(immdata$data, target, .col = "aa")
vis(tc, .plot = "smooth", .order=c(2,1,3,4,5))


query <- subset(df1, CDR3.aa %in% target) #subset df from mixcr with public clone list


##_-------------------------------CIRCOS------------------------------------##

file_list_int <- list.files(file_path)

file_list <- file_list_int[file_list_int != "metadata.txt"]

for (file in file_list){
  
  
  try <- file
  depth <- read.table(paste(file_path,"/",try,sep=""),header=T,sep="\t")
  
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
  
  dev.copy2pdf(file=paste(file_path,file,".pdf",sep=""))
  
}

##----------------------------Annotation----------------------------------------------###

tbadb = dbLoad("C:/Export/TCR_db/TBAdb.xlsx", "tbadb", .species = "Homo Sapiens", .chain = c("TRG"))

anno_q <- dbAnnotate(immdata$data, tbadb, "CDR3.aa", "CDR3.gamma.aa")

query_list <- anno_q[[1]]

anno_query <- subset(tbadb, CDR3.gamma.aa %in% query_list) 



##-----------------------------------------------------###


p <- ggplot(df, aes(x = timepoint,group=1))
p <- p + geom_line(aes(y = Copies, colour = "Copies/mg DNA"))

p <- p + geom_line(aes(y = Clones, colour = "Clones"))
p <- p + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Clone Abundance"))


p <- p + scale_x_discrete(limits = timepoint)

p <- p + scale_colour_manual(values = c("blue", "red"))+ theme_bw()

p <- p + theme(legend.position = c(0.3, 0.9)) + ggtitle("CATWDRPYYKKLF-TRGV3/TRGJ1 Abundance \nby CAR qPCR ")


#---------immunoSEQ Export Descriptive Metrics-------------#

datalist= list()


immunSEQsum <- function(file_path){
  
  df_list <- list.files(file_path)
  
  for (df in df_list){
  
  samp <- read.table(paste(file_path,"/",df,sep=""), header=TRUE, sep="\t")  ##### Necessary to convert clone ID count to true integer***
  
  template_num <- sum(samp$count..templates.reads.)
  
  query <- subset(samp, sequenceStatus=='In')
  
  query$prodFreq <- query$count..templates.reads./sum(query$count..templates.reads.)
  
  mxProd <- max(query$prodFreq) * 100 ## <- max prod  freq
  
  entro <- entropy(query$count..templates.reads.)
  
  dominance <- diversity(query$prodFreq, index = "simpson", MARGIN = 1, base = exp(1)) ### <-- simpson's 
  
  simP_clone <- round(sqrt(1-dominance),4) ## <- Simpson Clonality
  
  prodRE <- sum(query$count..templates.reads.)
  
  Entropy <- round(entro,3)
  Clonality <- round(simP_clone,3)
  Max_Frequency <- round(mxProd,3)
  Gene_Rearrangements <- template_num
  
  
  
  overview <<- data.frame(Entropy,Clonality,Max_Frequency,Gene_Rearrangements,prodRE)
  
  datalist[[df]] <- overview # add it to your list
  
  }
  
  big_data <<- do.call(rbind, datalist)
  
  write.table(big_data, file=paste(file_path,"/","immunoSeqSum.tsv",sep=""), sep="\t",row.names=TRUE,quote = FALSE,col.names=TRUE)
  
  return(big_data) # mxProd, entro, dominance, simP_clone)
  
  
}


###---------------IGH Reptoire Analysis to determine Onco-clonotypes and MRD---------------------------##

df_list <- c("base","tp")


datalist= list()

datalist= list()


Bcell_immunSEQsum <- function(file_path,df_list){
  for (df in df_list){
    
    base <- read.table(paste(file_path,"/",df_list[1],sep=""), header=TRUE, sep="\t")  ##### Necessary to convert clone ID count to true integer***
    
    base_sort <- base[order(-base$productive_frequency),]
    
    top_clone <- as.character(base_sort$amino_acid[1])
    
    genome_int <- as.numeric(base_sort$templates[1])
    
    leuk_clone <- nrow(subset(base_sort, frequency > 0.10))
    
    leuk_prod_clone <- nrow(subset(base_sort, frequency > 0.10 & frame_type=='In'))
    
    base_dj <-as.character(subset(base[order(-base$templates),], rearrangement_type =='DJ')$rearrangement[1])

    ##--------------------------##
    unsort_samp <- read.table(paste(file_path,"/",df,sep=""), header=TRUE, sep="\t")  ##### Necessary to convert clone ID count to true integer***
    print(df)
    print(as.character(unsort_samp$counting_method[1]))

    Sample_Type <- as.character(unsort_samp$sample_catalog_tags[1])
    
    #-----Top CLone Metrics----#
    

    clone_set <- subset(unsort_samp, amino_acid %in% top_clone) 
    

    if (nrow(clone_set)==0)
    {
      genome<-"ND"
      mass_frac<-"N/A"
      log_red <- "N/A"}else{
    
    clone_sort <- clone_set[order(-clone_set$productive_frequency),]
    
    genome <- as.numeric(clone_sort$templates[1])
    
    mass_frac <- as.numeric(clone_sort$templates[1]/clone_sort$sample_cells_mass_estimate[1])
    
    if (as.character(genome)==as.character(genome_int))
    {log_red <- "N/A"}else{log_red <- as.character(round(log10(genome_int/genome),2))}

      }
    
    #-----Metrics------#
    
    samp_sort <- unsort_samp[order(-unsort_samp$templates),]
    
    Total_Clones <-sum(unsort_samp$template)
    
    DJ_clone <- subset(samp_sort, rearrangement_type =='DJ')$templates[1]
    DJ_clonotype <- if (base_dj == as.character(subset(samp_sort, rearrangement_type =='DJ')$rearrangement[1])){"Yes"}else{"No"}

    samp <- subset(samp_sort, rearrangement_type!='DJ')
    
    template_num <- sum(samp$templates)
    
    query <- subset(samp_sort, frame_type=='In')
    
    Clonotypes <- nrow(query) #rearrangements
    
    query$prodFreq <- query$seq_reads/sum(query$seq_reads)
    
    mxProd <- max(query$prodFreq) * 100 ## <- max prod  freq
    
    entro <- entropy(query$seq_reads)
    
    dominance <- diversity(query$prodFreq, index = "simpson", MARGIN = 1, base = exp(1)) ### <-- simpson's 
    
    simP_clone <- round(sqrt(1-dominance),4) ## <- Simpson Clonality
    
    Productive_Clones <- sum(query$templates)
    
    Entropy <- round(entro,3)
    Clonality <- round(simP_clone,3)
    Max_Frequency <- round(mxProd,3)
    Gene_Clones <- template_num
    
    #Clone#
    CLL_Clonotype <- top_clone
    Fraction_of_Nucleated <- mass_frac
    Estimated_No_genomes <- genome
    Log_Reduction <- log_red
    
    overview <<- data.frame(Sample_Type,Entropy,Clonality,Max_Frequency,Total_Clones,Gene_Clones,Productive_Clones, 
                            Clonotypes,leuk_clone,leuk_prod_clone,DJ_clone,DJ_clonotype,CLL_Clonotype,Fraction_of_Nucleated,Estimated_No_genomes,Log_Reduction)
    
    datalist[[df]] <- overview # add it to your list
    
    
    
  }
  
  big_data <<- do.call(rbind, datalist)

  write.table(big_data, file=paste(file_path,"/","immunoSeqSum_18415-0.tsv",sep=""), sep="\t",row.names=TRUE,quote = FALSE,col.names=TRUE)

  return(big_data) # mxProd, entro, dominance, simP_clone)
  
  
}


