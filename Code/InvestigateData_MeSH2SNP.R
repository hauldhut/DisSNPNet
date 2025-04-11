#################
#MeSH-SNP associations
library(hash)
library(ggplot2)

df.Assoc = NULL
for(chr in 1:22){
  df.MeSHSNP <- read.csv(paste0("~/Data/GWAS/CAUSALdb/Chr_",chr,"_Assoc.txt_BinaryInteraction.csv"), header = TRUE)
  NuAssoc = nrow(df.MeSHSNP)
  NuMeSH = length(unique(df.MeSHSNP$disease))
  NuSNP = length(unique(df.MeSHSNP$enhancer))
  df.Assoc = rbind(df.Assoc,data.frame(chr=chr, NuMeSH=NuMeSH, NuSNP=NuSNP,NuAssoc=NuAssoc))
}
dim(df.Assoc)

result <- cor.test(df.Assoc$NuSNP, df.Assoc$NuAssoc, method = "pearson")  # For Pearson correlation
# You can also use method = "pearson" or "spearman" or "kendall" for other types
# To view the results
print(result)

fontsize = 12

FigureFile = "../Figure/Figure_MeSH2SNP"

write.csv(df.Assoc,paste0(FigureFile,".csv"),row.names = FALSE)

df.Assoc$chr = as.factor(df.Assoc$chr)
h.para2plot = hash()
para_label = ""
for(para in colnames(df.Assoc)){
  if(para == "NuMeSH"){
    df.Assoc$para = df.Assoc$NuMeSH
    para_label = "Number of MeSH terms"
  }else if(para == "NuSNP"){
    df.Assoc$para = df.Assoc$NuSNP
    para_label = "Number of SNPs"
  }else if(para == "NuAssoc"){
    df.Assoc$para = df.Assoc$NuAssoc    
    para_label = "Number of MeSH-SNP Associations"
  }else{
    next
  }
  
  
  p = ggplot(data=df.Assoc, aes(x=chr, y=para)) +
    geom_bar(stat="identity",position=position_dodge())
  
  p = p + labs(title=para, x="\nChr", y = paste0(para_label,"\n")) +
    theme_light() +
    scale_fill_brewer(palette="Reds") +
    scale_y_continuous(minor_breaks = seq(min(df.Assoc$para), max(df.Assoc$para), by = (max(df.Assoc$para)-min(df.Assoc$para))/10)) +
    theme(
      text = element_text(size=fontsize),# All font sizes
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size=fontsize),
      # legend.title=element_blank(),#Remove legend title (Network)
      axis.text = element_text(size = fontsize),
      axis.title = element_text(size = fontsize)
    )
  h.para2plot[[para]] = p
  saveRDS(p,paste0(FigureFile,"_",para,".rdata"))#readRDS
}
para = "NuMeSH"
p = h.para2plot[[para]]

# #Put all together
# df.AssocByPara = NULL
# for(para in colnames(df.Assoc)){
#   if(para == "NuMeSH"){
#     df.Assoc$para = df.Assoc$NuMeSH
#   }else if(para == "NuSNP"){
#     df.Assoc$para = df.Assoc$NuSNP
#   }else if(para == "NuAssoc"){
#     df.Assoc$para = df.Assoc$NuAssoc    
#   }else{
#     next
#   }
#   df.AssocByPara = rbind(df.AssocByPara, data.frame(chr = df.Assoc$chr, value = df.Assoc$para, group = rep(para,nrow(df.Assoc))))  
# }
# df.AssocByPara
# 
# df.AssocByPara$chr = as.factor(df.AssocByPara$chr)
# pAll = ggplot(data=df.AssocByPara, aes(x=chr, y=value, fill=group)) +
#   geom_bar(stat="identity",position=position_dodge())
# 
# pAll = pAll + labs(title=para, x="\nChr", y = paste0("NuMeSH/NuSNP/NuAssoc","\n")) +
#   scale_fill_brewer(palette="Blues") +
#   theme_light() +
#   theme(
#     text = element_text(size=fontsize),# All font sizes
#     plot.title = element_text(hjust = 0.5),
#     legend.text = element_text(size=fontsize),
#     # legend.title=element_blank(),#Remove legend title (Network)
#     axis.text = element_text(size = fontsize),
#     axis.title = element_text(size = fontsize)
#   )
# 
# saveRDS(pAll,paste0(FigureFile,"_",para,".rdata"))#readRDS

library('cowplot')
# plot_grid(pAll,h.para2plot[["NuMeSH"]], h.para2plot[["NuSNP"]],h.para2plot[["NuAssoc"]], labels=c("A", "B", "C","D"), ncol = 2, nrow = 2)

# plot_grid(h.para2plot[["NuMeSH"]], h.para2plot[["NuSNP"]],h.para2plot[["NuAssoc"]], labels=c("A", "B", "C","D"), ncol = 2, nrow = 2)
plot_grid(h.para2plot[["NuMeSH"]], h.para2plot[["NuSNP"]],h.para2plot[["NuAssoc"]], labels=c("A", "B", "C"), ncol = 3, nrow = 1)

ggsave(paste0(FigureFile,".pdf"), width = 15, height = 5)

#################
#MeSH-SNP associations by Pop
Pops = c("EUR","EAS","AFR","SAS","AMR")
DataFolder = "~/Data/GWAS/CAUSALdb/"
df.PopNuAssoc = NULL
for(pop in Pops){
  PopAssocFile = paste0(DataFolder,pop,"_Assoc.txt")
  print(PopAssocFile)
  df.PopAssoc = read.delim(PopAssocFile, header = FALSE, sep = "\t")
  NuMeSH = nrow(df.PopAssoc)
  
  NuAssoc = 0
  AllSNPvec = vector()
  for(i in 1:NuMeSH){
    SNPvec = strsplit(df.PopAssoc[i,3],", ")[[1]]
    AllSNPvec = append(AllSNPvec,SNPvec)
    NuAssoc = NuAssoc + length(SNPvec)
  }
  
  AllSNPvec = unique(AllSNPvec)
  
  df.PopNuAssoc = rbind(df.PopNuAssoc,data.frame(Pop = pop,NuMeSH =NuMeSH,NuSNP = length(AllSNPvec), NuAssoc=NuAssoc))
}
dim(df.PopNuAssoc)
head(df.PopNuAssoc)

bp<- ggplot(df.PopNuAssoc, aes(x="", y=NuAssoc, fill=Pop))+
  geom_bar(width = 1, stat = "identity")

bp = bp + coord_polar("y", start=0)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    legend.position = c(1.1, 0.5),# Position legend inside the plot (closer to center)
    legend.justification = c("right", "center"),  # Anchor the legend's right-center to the position
    plot.margin = unit(c(0.5, 2, 1, 1), "lines")  # Add right margin to make space for the legend
  )

bp = bp + blank_theme +
  theme(axis.text.x=element_blank())
# geom_text(aes(label = NuAssoc))
bp

library("ggpmisc")
bp = bp +                                               # Add table to ggplot2 plot
  annotate(geom = "table",
           x = 0.5,
           y = -0.5,
           label = list(df.PopNuAssoc))

ggsave("../Figure/Figure_MeSH2SNP_Pop.pdf", width = 5, height = 5)

plot_grid(bp,h.para2plot[["NuMeSH"]], h.para2plot[["NuSNP"]],h.para2plot[["NuAssoc"]], labels=c("A", "B", "C","D"), ncol = 2, nrow = 2)

ggsave(paste0(FigureFile,"_New.pdf"), width = 10, height = 10)

################
#MeSHID_Net.txt
disease <- read.delim("../Data/MeSHID_Net.txt",header = FALSE)
head(disease)
disease.frame <- data.frame(disease[[1]], disease[[3]])
disease.weight = disease[[2]]

disease.g <- graph.data.frame(d = disease.frame, directed = FALSE)
#add weight of the disease net
E(disease.g)$weight <- disease.weight








