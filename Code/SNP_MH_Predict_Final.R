#Load all EFO 2 DOID mappings from EFO.obo
library('ontologyIndex')
library(hash)

# efo = get_OBO(file = "~/Data/Ontology/efo.obo", extract_tags = "everything")
efo = get_ontology(file = "~/Data/Ontology/efo.obo", extract_tags = "everything")

EFOIDVec = vector()
h.efoid2meshid = hash()
for(i in 1:length(efo$id)){
  id = efo$id[[i]]
  if(substr(id,1,3)=="EFO"){
    EFOIDVec = append(EFOIDVec,id)
    
    xrefvec = efo$xref[[id]]
    if(length(xrefvec)>0){
      meshidvec = vector()
      for(xref in xrefvec){
        pos = which(strsplit(xref, "")[[1]]==":")[1]
        if(is.na(pos)){ 
          print(xref)
          next
        }
        if(substr(xref,1,pos-1)=="MeSH" || substr(xref,1,pos-1)=="MESH" || substr(xref,1,pos-1)=="MSH"){
          meshid = substr(xref,pos+1,nchar(xref))
          meshidvec = append(meshidvec,meshid)
        }
      }
      if(length(meshidvec)>0){
        h.efoid2meshid[[id]] = meshidvec  
      }
    }
  }
}
EFOIDVec = unique(EFOIDVec)
length(EFOIDVec)
length(h.efoid2meshid)

#Functions
removeSmallNetCluster <- function(g, size){
  dg <- decompose.graph(g) 
  
  RemovedNodeList = vector()
  for(ci in 1:length(dg)){
    nodes = V(dg[[ci]])
    if(length(nodes)<=size){
      for(ni in 1:length(nodes)){
        RemovedNodeList = append(RemovedNodeList, nodes[[ni]]$name)  
      }
    }
  }
  
  RemovedNodeList = unique(RemovedNodeList)
  g2 = delete_vertices(g, RemovedNodeList)
  
  return(g2)
}

Method = "H"

start_time <- Sys.time()

library(RandomWalkRestartMH)
library(igraph)
library(foreach)
library(doParallel)
library(ROCR)

setwd("~/Manuscripts/70MeSHSNPAssoc/Code")

phase = 3
ld = "r2def.ld"#r2def.ld/r208.ld

df.NetProp = NULL
for(chr in 1:22){
  snpnetfile = paste0("~/Data/1KGP/Phase",phase,"/1kg_phase",phase,"_chr",chr,"_",ld,"_Net.txt")
  net <- read.table(snpnetfile, header = FALSE) 
  netnode = union(net[[1]],net[[3]])
  df.NetProp = rbind(df.NetProp, data.frame(chr = chr, NuNetNode = length(netnode)))
}

df.NetProp = df.NetProp[order(df.NetProp$NuNetNode),]
df.NetProp

size = 2
# for(chri in 1:nrow(df.NetProp)){#nrow(df.NetProp)
for(chr in c(10,2)){  
  # chri = 1
  
  # chr=df.NetProp$chr[chri]
  # cat("PROCESSING: chri=",chri,"chr=",chr,"\n")
  
  snpnetfile = paste0("~/Data/1KGP/Phase",phase,"/1kg_phase",phase,"_chr",chr,"_",ld, "_Net.txt")#r2 = 0.8
  
  print(snpnetfile)
  enhancer1 <- read.delim(snpnetfile,header = FALSE)
  
  enhancer1.frame <- data.frame(enhancer1[[1]], enhancer1[[3]])
  enhancer1.g <- graph.data.frame(d = enhancer1.frame, directed = FALSE)
  enhancer1.weight = enhancer1[[2]]
  E(enhancer1.g)$weight <- enhancer1.weight
  
  print("Before")
  print(length(V(enhancer1.g)))
  print(length(E(enhancer1.g)))
  
  enhancer1.g = removeSmallNetCluster(enhancer1.g, size = size)
  
  print("After")
  print(length(V(enhancer1.g)))
  print(length(E(enhancer1.g)))
  
  enhancer_MultiplexObject <- create.multiplex(list(enhancer1.g),Layers_Name = c("enhancer"))
  
  AdjMatrix_enhancer <- compute.adjacency.matrix(enhancer_MultiplexObject)
  AdjMatrixNorm_enhancer <- normalize.multiplex.adjacency(AdjMatrix_enhancer)
  
  #Add disease nw
  #MeSHID_Net.txt
  disease <- read.delim("../Data/MeSHID_Net.txt",header = FALSE)
  head(disease)
  disease.frame <- data.frame(disease[[1]], disease[[3]])
  disease.weight = disease[[2]]
  
  disease.g <- graph.data.frame(d = disease.frame, directed = FALSE)
  #add weight of the disease net
  E(disease.g)$weight <- disease.weight
  
  disease_MultiplexObject <- create.multiplex(list(disease.g),
                                              Layers_Name = c("disease"))
  
  #Add EDRelation
  
  ED.frame <- read.csv(paste0("~/Data/GWAS/CAUSALdb/Chr_",chr,"_Assoc.txt_BinaryInteraction.csv"), header = TRUE)
  
  ED.frame <- ED.frame[which(ED.frame$enhancer %in% enhancer_MultiplexObject$Pool_of_Nodes),]
  ED.frame <- ED.frame[which(ED.frame$disease %in% disease_MultiplexObject$Pool_of_Nodes),]
  head(ED.frame)
  dim(ED.frame)
  
  #Create multiplex-heterosgenous nw
  enhancer_Disease_Net <- create.multiplexHet(enhancer_MultiplexObject, disease_MultiplexObject, 
                                              ED.frame)
  
  enhancerHetTranMatrix <- compute.transition.matrix(enhancer_Disease_Net)
  
  #pick unique disease
  unique_disease <- unique(ED.frame$disease)
  
  #Extract DOID2TraitMap
  
  MeSHInfoFile = "~/Data/GWAS/CAUSALdb/MeSHInfo.txt"
  MeSHInfo = read.delim(MeSHInfoFile, sep = "\t", header = FALSE)
  
  MeSHID2TermMap = hash()
  MeSHTerm2IDMap = hash()
  MeSHID2TraitMap = hash()
  
  for(i in 1:nrow(MeSHInfo)){
    meshid = MeSHInfo[i,1]
    meshterm = MeSHInfo[i,2]
    MeSHID2TermMap[[meshid]] = meshterm
    MeSHTerm2IDMap[[meshterm]] = meshid
    traitvec = strsplit(MeSHInfo[i,3],"//")[[1]]
    MeSHID2TraitMap[[meshid]] = traitvec
  }
  length(MeSHID2TermMap)
  length(MeSHTerm2IDMap)
  length(MeSHID2TraitMap)
  
  #loop through
  
  h.d2ranking = hash()
  
  for(di in 1:length(unique_disease)){ 
    # di = 1

    
    seeddisease = unique_disease[[di]]
    
    cat("==> Disease ", di, "/", length(unique_disease), ": ", seeddisease,"\n")
    
    disease_relation = ED.frame[which(ED.frame$disease==seeddisease),]
    
    SeedEnhancer = disease_relation$enhancer[c(which(disease_relation$disease==seeddisease))]
    
    #compute 
    RWRH_enhancer_Disease_Results <- Random.Walk.Restart.MultiplexHet(enhancerHetTranMatrix,
                                                                      enhancer_Disease_Net,SeedEnhancer,
                                                                      seeddisease, r = 0.5)
    
    tf = RWRH_enhancer_Disease_Results$RWRMH_Multiplex1
    h.d2ranking[[seeddisease]] = tf  
  }
  
  #####################################################
  
  #Summarize for each k
  maxk=30
  h.K2TotalEvidencedNewAssoc = hash()
  h.K2topKEnhEvidence = hash()
  
  #Store top k predictions
  for(k in seq(10, maxk, by = 10)){
    
    cat("k =",k,"\n")
    df.topKEnh <- data.frame(matrix(ncol = 2+k, nrow = 0))
    
    for(meshid in unique_disease){ 
      tf = h.d2ranking[[meshid]]
      df.topKEnh[nrow(df.topKEnh)+1,] = c(as.character(meshid),as.character(MeSHID2TermMap[[meshid]]),as.vector(tf$NodeNames[1:k]))
    }
    topKfile = paste0("../Results/Prediction/",Method,"_byDisease_Phase",phase,"_Chr",chr,"_",ld,"_predict_top",k,".txt")
    write.table(df.topKEnh, topKfile, na ="", row.names=FALSE, col.names = FALSE, sep='\t', quote=FALSE)
  }
    
  #Find evidence for top k
  df.topKEnhEvidence = NULL
  
  for(k in seq(10, maxk, by = 10)){
    cat("SEARCH FOR K in range:",(k-10)+1,":",k,"\n")
    library('phenoscanner')
    
    h.meshid2evid = hash()
    di=0
    
    
    topKEvidfile = paste0("../Results/Prediction/",Method,"_byDisease_Phase",phase,"_Chr",chr,"_",ld,"_predict_top",k,"_evid.txt")
    if(file.exists(topKEvidfile)){
      
      df.topKEnhEvidence = read.delim(topKEvidfile, header = FALSE, sep = "\t")
      cat("DATA FOR",k, "EXISTED!",dim(df.topKEnhEvidence),"\n")
      next
    }
    df.topKEnhEvidenceByKrange = NULL  #1:10/11:20,21:30,...,(k-10)+1:k
    for(meshid in unique_disease){ 
      
      di = di+1
      # if(di<94) next
      cat(meshid,"\n")
      
      meshterm = MeSHID2TermMap[[meshid]]
      meshtraitvec = MeSHID2TraitMap[[meshid]]
      
      tf = h.d2ranking[[meshid]]
      
      starti = (k-10)+1
      snpquery = as.vector(tf$NodeNames[starti:k])
      # print(length(snpquery))
      res <- phenoscanner(snpquery=snpquery)
      
      EnhInfo = res$results
      EnhInfo = EnhInfo[EnhInfo$rsid %in% snpquery,]
      
      if(nrow(EnhInfo)==0) next
        
      df.topKEnhEvidenceByDisease = NULL
      for(i in 1:nrow(EnhInfo)){
        trait = EnhInfo$trait[i]
        # cat(trait,"\n")
        efoid = EnhInfo$efo[i]
        efoid = gsub("_", ":",efoid)
        
        GWASMeSHIDSet = h.efoid2meshid[[efoid]]
        if(meshid %in% GWASMeSHIDSet){
          # cat(meshid,"\t",meshterm,"\t",efoid,"\t",trait,"\t",EnhInfo$rsid[i],"(",EnhInfo$p[i],")","\t",EnhInfo$pmid[i],"\n")
          df.topKEnhEvidenceByDisease = rbind(df.topKEnhEvidenceByDisease, data.frame(meshid = meshid,meshterm=meshterm, trait=trait, efoid=efoid,rsid = EnhInfo$rsid[i],p=EnhInfo$p[i], pmid=EnhInfo$pmid[i]))
        }
      }
      
      nevid = 0
      if(!is.null(df.topKEnhEvidenceByDisease)){
        # h.meshid2evid[[meshid]] = df.topKEnhEvidenceByDisease  
        nevid = nrow(df.topKEnhEvidenceByDisease)
        df.topKEnhEvidenceByKrange = rbind(df.topKEnhEvidenceByKrange, df.topKEnhEvidenceByDisease)
      }
      cat("SEARCH FOR K in range:",(k-10)+1,":",k,"\n")
      cat(chr,di,"/", length(unique_disease),":",meshid,nrow(EnhInfo),nevid,"\n")
    }#end for(meshid in unique_disease){ 
    # length(h.meshid2evid)
    
    if(!is.null(df.topKEnhEvidence)){
      colnames(df.topKEnhEvidence)=colnames(df.topKEnhEvidenceByKrange)  
    }
    df.topKEnhEvidence = rbind(df.topKEnhEvidence, df.topKEnhEvidenceByKrange)
    
    dim(df.topKEnhEvidence)
    
    write.table(df.topKEnhEvidence, topKEvidfile, na ="", row.names=FALSE, col.names = FALSE, sep='\t', quote=FALSE)
    
  }
  print("===============")
  print("Waiting for 60 mins...!")
  Sys.sleep(60*60)
}

# #Test
# meshid = "D050723"#D002769/D003069
# meshterm = MeSHID2TermMap[[meshid]]
# meshtraitvec = MeSHID2TraitMap[[meshid]]
# 
# tf = h.d2ranking[[meshid]]
# 
# snpquery = as.vector(tf$NodeNames[1:k])
# res <- phenoscanner(snpquery=snpquery)
# 
# EnhInfo = res$results
# dim(EnhInfo)
# 
# df.topKEnhEvidence = NULL
# for(i in 1:nrow(EnhInfo)){
#   trait = EnhInfo$trait[i]
#   
#   efoid = EnhInfo$efo[i]
#   efoid = gsub("_", ":",efoid)
#   
#   cat(efoid,"\n")
#   
#   GWASMeSHIDSet = h.efoid2meshid[[efoid]]
#   if(meshid %in% GWASMeSHIDSet){
#     cat(meshid,"\t",meshterm,"\t",efoid,"\t",trait,"\t",EnhInfo$rsid[i],"(",EnhInfo$p[i],")","\t",EnhInfo$pmid[i],"\n")
#     # df.topKEnhEvidence = rbind(df.topKEnhEvidence, data.frame(meshid = meshid,meshterm=meshterm, trait=trait, efoid=efoid,rsid = paste0(EnhInfo$rsid[i]," (",EnhInfo$p[i],")"), pmid=EnhInfo$pmid[i]))
#   }
#   
# }
