library(hash)
library(ggplot2)
library(RColorBrewer)
library(ROCR)

setwd("~/Manuscripts/70MeSHSNPAssoc/Code")

Method = "H"#M/H

Methods =c("M","H")

Phases = c(1,3)
LDs = c("r208","r2def")

h.ld2Thres = hash()
h.ld2Thres[["r208"]]="0.8"
h.ld2Thres[["r2def"]]="0.2"

SummaryType = paste0("../Figure/Summary_topKEvidence_",Method)

h.config2topKEvid = hash()

phase = 1
ld = "r2def"

for(phase in Phases){
  for(ld in LDs){
    
    config = paste0("Phase",phase,"_",ld)
    
    cat("PROCESSING FOR", config, "\n")
    
    #Load New Results (traverse from chri 1 -> 22 (chr ordered by NuNode increasingly), increase size (size of clusters to eliminate) one whenever an error occurs)
    df.topKEvid <-NULL
    for(chr in 1:22){
      cat("==>chr", chr, "\n")
      
      maxk = 30
      h.k2NuEvid = hash()
      for(k in seq(10, maxk, by = 10)){
        topKEvidfile = paste0("../Results/Prediction/",Method,"_byDisease_Phase",phase,"_Chr",chr,"_", ld,".ld_predict_top",k,"_evid.txt")
        print(topKEvidfile)
        
        topKEvid = NaN
        if(file.exists(topKEvidfile)){
          res <- read.delim(topKEvidfile, header = FALSE, sep = "\t") 
          topKEvid = nrow(res)
          
        }else{
          print("File is not existing...!")
        }
        h.k2NuEvid[[as.character(k)]] = topKEvid
      }
      cat("chr =", chr, "top10 =", h.k2NuEvid[["10"]], "top20 =", h.k2NuEvid[["20"]], "top30 =",h.k2NuEvid[["30"]],"\n")
      df.topKEvid <- rbind(df.topKEvid, data.frame(chr = chr, top10 = h.k2NuEvid[["10"]], top20 = h.k2NuEvid[["20"]], top30=h.k2NuEvid[["30"]]))
    }
    
    SummaryFile = paste0(SummaryType,"_",config,".csv")
    write.csv(df.topKEvid,SummaryFile, row.names = FALSE, quote = FALSE)
    
    h.config2topKEvid[[config]] = df.topKEvid
  }
}
length(h.config2topKEvid)
h.config2topKEvid[["Phase1_r208"]]

h.config2topKEvid_temp = h.config2topKEvid

#Load summarized AUROC/AUPR for each config
SummaryType = paste0("../Figure/Summary_topKEvidence_",Method)

h.config2topKEvid = hash()
for(phase in Phases){
  for(ld in LDs){
    config = paste0("Phase",phase,"_",ld)
    SummaryFile = paste0(SummaryType,"_",config,".csv")
    df.topKEvid = read.csv(SummaryFile, header = TRUE)
    h.config2topKEvid[[config]]=df.topKEvid
  }
}
length(h.config2topKEvid)

#Draw
h.config2topKEvidPlot = hash()
for(phase in Phases){
  for(ld in LDs){
    config = paste0("Phase",phase,"_",ld)
    df.topKEvid = h.config2topKEvid[[config]]
    
    df.topK = data.frame(chr=df.topKEvid$chr, NuEvid = df.topKEvid$top10, topK = rep("top10",nrow(df.topKEvid)))
    df.topK = rbind(df.topK, data.frame(chr=df.topKEvid$chr, NuEvid = df.topKEvid$top20, topK = rep("top20",nrow(df.topKEvid))))
    df.topK = rbind(df.topK, data.frame(chr=df.topKEvid$chr, NuEvid = df.topKEvid$top30, topK = rep("top30",nrow(df.topKEvid))))

    cat(config, dim(df.topK),"\n")
    
    df.topK$chr = as.factor(df.topK$chr)

    fontsize = 16
    p = ggplot(data=df.topK, aes(x=chr, y=NuEvid, fill=topK)) +
      geom_bar(stat="identity",position=position_dodge())

    p = p + labs(title=paste0("Phase=", phase, ", LD=", h.ld2Thres[[ld]]), x="\nChr", y = "NuEvid") +
      # scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
      theme_light() +
      scale_fill_brewer(palette="Reds") +
      theme(
        panel.background = element_rect(fill = 'white', colour = 'white'),
        text = element_text(size=fontsize),# All font sizes
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=fontsize),
        # legend.title=element_blank(),#Remove legend title (Network)
        axis.text = element_text(size = fontsize),
        axis.title = element_text(size = fontsize)
      )

    h.config2topKEvidPlot[[config]] = p
  }
}
length(h.config2topKEvidPlot)
# h.config2topKEvidPlot[["Phase1_r2def.ld"]]
# h.config2topKEvidPlot[["Phase1_r208.ld"]]
# h.config2topKEvidPlot[["Phase3_r2def.ld"]]
# h.config2topKEvidPlot[["Phase3_r208.ld"]]

library('cowplot')
# pA = h.config2topKEvidPlot[["Phase1_r2def"]] + theme(legend.position = "none",plot.margin = unit(c(1, 3, 1, 1), "lines"))
# pB = h.config2topKEvidPlot[["Phase1_r208"]] + labs(y = NULL)
# pC = h.config2topKEvidPlot[["Phase3_r2def"]] + theme(legend.position = "none",plot.margin = unit(c(1, 3, 1, 1), "lines"))
# pD = h.config2topKEvidPlot[["Phase3_r208"]] + labs(y = NULL)
# 
# plot_grid(pA, pB, pC, pD, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2)
# 
# FigureFile = paste0("../Figure/Summary_topKEvidence_",Method,"_New_v3.pdf")
# ggsave(FigureFile, width = 15, height = 10)


# pB = h.config2topKEvidPlot[["Phase1_r208"]] + theme(legend.position = "none",plot.margin = unit(c(1, 3, 1, 1), "lines"))
# pD = h.config2topKEvidPlot[["Phase3_r208"]] + labs(y = NULL)
# 
# plot_grid(pB, pD, labels=c("A", "B"), ncol = 2, nrow = 1)
# 
# FigureFile = paste0("../Figure/Summary_topKEvidence_",Method,"_New_v3.pdf")
# ggsave(FigureFile, width = 15, height = 5)


pD = h.config2topKEvidPlot[["Phase3_r208"]]

plot_grid(pD, ncol = 1, nrow = 1)

FigureFile = paste0("../Figure/Summary_topKEvidence_",Method,"_New_v3.pdf")
ggsave(FigureFile, width = 15, height = 5)

