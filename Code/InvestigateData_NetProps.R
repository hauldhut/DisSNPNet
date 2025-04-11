  library(hash)
  library(ggplot2)
  library(RColorBrewer)
  library("patchwork")
  library(igraph)
  phase = 3
  ld = "r2def" 
  
  setwd("~/Manuscripts/70MeSHSNPAssoc/Code")
  
  FigureFile = "../Figure/Figure_NetProp_"
  ################
  #LD Nets
  Phases = c(1,3)
  LDs = c("r208","r2def")
  h.ld2Thres = hash()
  h.ld2Thres[["r208"]]="0.8"
  h.ld2Thres[["r2def"]]="0.2"
  
  h.config2dfNetProp = hash()
  for(phase in Phases){
    for(ld in LDs){
      config = paste0("Phase",phase,"_",ld)
      df.NetProp = NULL
      for(chr in 1:22){
        snpnetfile = paste0("~/Data/1KGP/Phase",phase,"/1kg_phase",phase,"_chr",chr,"_",ld,".ld_Net.txt")
        print(snpnetfile)
        net <- read.table(snpnetfile, header = FALSE) 
        net.frame <- data.frame(net[[1]], net[[3]])
        net.g <- graph.data.frame(d = net.frame, directed = FALSE)
        net.weight = net[[2]]
        E(net.g)$weight <- net.weight
        
        cl = clusters(net.g)
        minClSize = min(cl$csize)
        maxClSize = max(cl$csize)
        nuCl = cl$no
        
        dim(net)
        netnode = union(net[[1]],net[[3]])
        
        NuNetLink = nrow(net)
        NuNetNode = length(netnode)
        
        cat("==>",phase,ld,chr,NuNetNode,nuCl, NuNetLink, minClSize, maxClSize,"\n")
        df.NetProp <- rbind(df.NetProp,data.frame(chr=chr,NuNetNode=NuNetNode, NuNetLink=NuNetLink,nuCl=nuCl, minClSize = minClSize, maxClSize = maxClSize))  
      }
      h.config2dfNetProp[[config]] = df.NetProp
    }
  }
  length(h.config2dfNetProp)
  # h.config2dfNetProp[["Phase1_r2def"]]
  # h.config2dfNetProp[["Phase1_r208"]]
  # h.config2dfNetProp[["Phase3_r2def"]]
  # h.config2dfNetProp[["Phase4_r208"]]
  
  for(phase in Phases){
    for(ld in LDs){
      config = paste0("Phase",phase,"_",ld)
      write.csv(h.config2dfNetProp[[config]][,-5],paste0(FigureFile,config,".csv"),row.names = FALSE)
    }
  }
  
  h.paraconfig2plot = hash()
  
  
  
  for(phase in Phases){
    for(ld in LDs){
      config = paste0("Phase",phase,"_",ld)
      cat("PROCESSING FOR ", config, "\n")
      
      df.NetProp = h.config2dfNetProp[[config]]
      
      
      df.NetProp$chr = as.factor(df.NetProp$chr)
      # h.para2plot = hash()
      for(para in colnames(df.NetProp)){
        
        para_label = ""
        if(para == "NuNetNode"){
          df.NetProp$para = df.NetProp$NuNetNode
          para_label = "Number of SNPs"
        }else if(para == "NuNetLink"){
          df.NetProp$para = df.NetProp$NuNetLink
          para_label = "Number of Links"
        }else if(para == "nuCl"){
          df.NetProp$para = df.NetProp$nuCl    
          para_label = "Number of Clusters"
        }else if(para == "maxClSize"){
          df.NetProp$para = df.NetProp$maxClSize
          para_label = "Maximum Cluster Size"
        }else{
          next
        }
        
        fontsize = 16
        
        p = ggplot(data=df.NetProp, aes(x=chr, y=para)) +
          geom_bar(stat="identity",position=position_dodge())
        
        p = p + labs(title=paste0("Phase=",phase,", LD=",h.ld2Thres[[ld]]), x="\nChr", y = paste0(para_label,"\n")) +
          theme_light() +
          scale_fill_brewer(palette="Reds") +
          scale_y_continuous(minor_breaks = seq(min(df.NetProp$para), max(df.NetProp$para), by = (max(df.NetProp$para)-min(df.NetProp$para))/10)) +
          theme(
            text = element_text(size=fontsize),# All font sizes
            plot.title = element_text(hjust = 0.5),
            legend.text = element_text(size=fontsize),
            # legend.title=element_blank(),#Remove legend title (Network)
            axis.text = element_text(size = fontsize),
            axis.title = element_text(size = fontsize)
          )
        # h.para2plot[[para]] = p
        paraconfig = paste0(para,"_",config)
        h.paraconfig2plot[[paraconfig]] = p
        saveRDS(p,paste0(FigureFile,paraconfig,".rdata"))#readRDS
        ggsave(paste0(FigureFile,paraconfig,".pdf"), width = 10, height = 5)
      }
      
    }
  }
  
  length(h.paraconfig2plot)
  
  library('cowplot')
  # plot_grid(h.para2plot[["NuNetNode"]], h.para2plot[["NuNetLink"]],h.para2plot[["nuCl"]],h.para2plot[["maxClSize"]], labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2)
  # 
  # ggsave(paste0(FigureFile,"_",config,"_New.pdf"), width = 10, height = 10)
  
  # para = "NuNetNode"
  # p1 = h.paraconfig2plot[[paste0(para,"_Phase1_r2def")]] + theme(legend.position = "none")
  # p2 = h.paraconfig2plot[[paste0(para,"_Phase1_r208")]] + theme(legend.position = "none")+ labs(y = "")
  # p3 = h.paraconfig2plot[[paste0(para,"_Phase3_r2def")]]+ theme(legend.position = "none")+ labs(y = "")
  # p4 = h.paraconfig2plot[[paste0(para,"_Phase3_r208")]] + labs(y = NULL)
  # 
  # para = "NuNetLink"
  # p5 = h.paraconfig2plot[[paste0(para,"_Phase1_r2def")]] + theme(legend.position = "none")
  # p6 = h.paraconfig2plot[[paste0(para,"_Phase1_r208")]] + theme(legend.position = "none")+ labs(y = "")
  # p7 = h.paraconfig2plot[[paste0(para,"_Phase3_r2def")]]+ theme(legend.position = "none")+ labs(y = "")
  # p8 = h.paraconfig2plot[[paste0(para,"_Phase3_r208")]] + labs(y = NULL)
  # 
  # para = "nuCl"
  # p9 = h.paraconfig2plot[[paste0(para,"_Phase1_r2def")]] + theme(legend.position = "none")
  # p10 = h.paraconfig2plot[[paste0(para,"_Phase1_r208")]] + theme(legend.position = "none")+ labs(y = "")
  # p11 = h.paraconfig2plot[[paste0(para,"_Phase3_r2def")]]+ theme(legend.position = "none")+ labs(y = "")
  # p12 = h.paraconfig2plot[[paste0(para,"_Phase3_r208")]] + labs(y = NULL)
  # 
  # para = "maxClSize"
  # p13 = h.paraconfig2plot[[paste0(para,"_Phase1_r2def")]] + theme(legend.position = "none")
  # p14 = h.paraconfig2plot[[paste0(para,"_Phase1_r208")]] + theme(legend.position = "none")+ labs(y = "")
  # p15 = h.paraconfig2plot[[paste0(para,"_Phase3_r2def")]]+ theme(legend.position = "none")+ labs(y = "")
  # p16 = h.paraconfig2plot[[paste0(para,"_Phase3_r208")]] + labs(y = NULL)
  # 
  # plot_grid(p1, p2, p3, p4, p5, p6, p7, p8,p9,p10,p11,p12,p13,p14,p15,p16, labels=c("A", "", "", "", "B", "", "", "", "C", "", "", "", "D", "", "", ""), ncol = 4, nrow = 4)
  # ggsave("../Figure/Figure_NetProps.pdf", width = 25, height = 20)
  
  para = "NuNetNode"
  p1 = h.paraconfig2plot[[paste0(para,"_Phase1_r2def")]] + theme(legend.position = "none",plot.margin = unit(c(1, 1, 1, 1), "lines")) 
  p2 = h.paraconfig2plot[[paste0(para,"_Phase1_r208")]] + labs(y = NULL)
  p3 = h.paraconfig2plot[[paste0(para,"_Phase3_r2def")]] + theme(legend.position = "none",plot.margin = unit(c(1, 1, 1, 1), "lines")) 
  p4 = h.paraconfig2plot[[paste0(para,"_Phase3_r208")]] + labs(y = NULL)
  
  para = "NuNetLink"
  p5 = h.paraconfig2plot[[paste0(para,"_Phase1_r2def")]] + theme(legend.position = "none",plot.margin = unit(c(1, 1, 1, 1), "lines")) 
  p6 = h.paraconfig2plot[[paste0(para,"_Phase1_r208")]] + labs(y = NULL)
  p7 = h.paraconfig2plot[[paste0(para,"_Phase3_r2def")]] + theme(legend.position = "none",plot.margin = unit(c(1, 1, 1, 1), "lines")) 
  p8 = h.paraconfig2plot[[paste0(para,"_Phase3_r208")]] + labs(y = NULL)
  
  para = "nuCl"
  p9 = h.paraconfig2plot[[paste0(para,"_Phase1_r2def")]] + theme(legend.position = "none",plot.margin = unit(c(1, 1, 1, 1), "lines")) 
  p10 = h.paraconfig2plot[[paste0(para,"_Phase1_r208")]] + labs(y = NULL)
  p11 = h.paraconfig2plot[[paste0(para,"_Phase3_r2def")]] + theme(legend.position = "none",plot.margin = unit(c(1, 1, 1, 1), "lines")) 
  p12 = h.paraconfig2plot[[paste0(para,"_Phase3_r208")]] + labs(y = NULL)
  
  para = "maxClSize"
  p13 = h.paraconfig2plot[[paste0(para,"_Phase1_r2def")]] + theme(legend.position = "none",plot.margin = unit(c(1, 1, 1, 1), "lines")) 
  p14 = h.paraconfig2plot[[paste0(para,"_Phase1_r208")]] + labs(y = NULL)
  p15 = h.paraconfig2plot[[paste0(para,"_Phase3_r2def")]] + theme(legend.position = "none",plot.margin = unit(c(1, 1, 1, 1), "lines")) 
  p16 = h.paraconfig2plot[[paste0(para,"_Phase3_r208")]] + labs(y = NULL)
  
  plot_grid(p1, p2, p5, p6, p9,p10,p13,p14, labels=c("A", "B", "C", "D", "E", "F", "G", "H"), ncol = 2, nrow = 4)
  ggsave("../Figure/Figure_NetProps_Phase1.pdf", width = 15, height = 20)
  plot_grid(p3, p4, p7, p8,p11,p12,p15,p16, labels=c("A", "B", "C", "D", "E", "F", "G", "H"), ncol = 2, nrow = 4)
  ggsave("../Figure/Figure_NetProps_Phase3.pdf", width = 15, height = 20)
