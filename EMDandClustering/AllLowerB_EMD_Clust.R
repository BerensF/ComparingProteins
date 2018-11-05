# ---------------------------------
# Created on: 11.08.2018
#     Author: Felix Berens
# 
# Functions that handle the calculation of all pairwise comparisons
# of the lower bounds (AllvsAll.LowerBound), and earth-mover's distances (AllvsAll.EMD).
# Also functions for the calculation of the isosurfaces (CreatALLdx) is here
# and the clustering (AllvsAll.Cluster).
# ---------------------------------

is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
}

if(!is.installed("emdist")){install.packages("emdist")}
if(!is.installed("ggplot2")){install.packages("ggplot2")}
if(!is.installed("ggdendro")){install.packages("ggdendro")}
if(!is.installed("plot3D")){install.packages("plot3D")}
if(!is.installed("emdist")){install.packages("emdist")}
if(!is.installed("tcltk")){install.packages("tcltk")}
library("cluster")
library("tcltk")
library("ggplot2")
library("ggdendro")
library("plot3D")
library("emdist")

rm(is.installed)

CreatALLdx <- function(ListOfProtNames,PathToProtData)
{
  for(i in 1:NROW(ListOfProtNames)){
    # needs only name without extension
    fileName <- ListOfProtNames[i]
    
    # needs full path
    inPath=paste(PathToProtData,"/",ListOfProtNames[i],sep="")
    
    # needs full path
    outPath=inPath
  
    if(!file.exists(paste(outPath,"/",fileName,"_pot_positive.pts",sep="")) | 
       !file.exists(paste(outPath,"/",fileName,"_pot_negative.pts",sep="")))
    {
      dxData <- read.csv(file=paste(inPath,"/",fileName,"_pot",".dx",sep=""),
                         sep=' ', skip= 11, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)
      dxData <- head(dxData,-10)
      
      v1 = as.numeric(dxData$V1)
      v2 = as.numeric(dxData$V2)
      v3 = as.numeric(dxData$V3)
      
      size = 129
      x <- c(1:size)
      y <- c(1:size)
      z <- c(1:size)
      
      merged <- as.vector(rbind(v1,v2,v3)) 
      V <- array(merged, c(size,size,size))
      
      # print("creating isosurface ...")
      iso <- createisosurf(x, y, z, V, level = 1.0)
      iso2 <- createisosurf(x, y, z, V, level = -1.0)
      
      iso <- iso[!duplicated(iso), ]
      iso2 <- iso2[!duplicated(iso2), ]
      # print(paste("writing to file ", outPath, "/", fileName, ".pts ...", sep = ""))
      write.table(iso, file = paste(outPath,"/",fileName,"_pot_positive.pts",sep=""),
                  row.names = F,na = "",sep = ";",dec = ".")
      write.table(iso2, file = paste(outPath,"/",fileName,"_pot_negative.pts",sep=""),
                  row.names = F,na = "",sep = ";",dec = ".")
   }
  }
}


AllvsAll.LowerBound <- function(ListOfProtNames,PathToProtData,PathToOutput,PathToProgram,n,m,pot,pb)
{
  total <- NROW(ListOfProtNames)^2/2+NROW(ListOfProtNames)/2
  a <- pot == "positive"
  count <- 0
  
  for(i in 1:NROW(ListOfProtNames))
  {
    OutFilePath <- paste(PathToOutput,"/",ListOfProtNames[i],"_",pot,"_",n,sep="")
    if(!dir.exists(OutFilePath))
    {
      dir.create(OutFilePath,recursive = T)
    }
    ProtAPath <- paste(PathToProtData,"/",ListOfProtNames[i],"/",ListOfProtNames[i],"_pot_",pot,".pts",sep="")
    for(j in i:NROW(ListOfProtNames))
    {
      
      ProtBPath <- paste(PathToProtData,"/",ListOfProtNames[j],"/",ListOfProtNames[j],"_pot_",pot,".pts",sep="")
      if(!file.exists(paste(OutFilePath,"/",ListOfProtNames[i],"_",ListOfProtNames[j],"_",pot,"_",n,sep="")))
      {
        time <- Sys.time()
        write.csv(as.numeric(system2(PathToProgram,c(ProtAPath,ProtBPath,n,m),stdout = T)),
                  paste(OutFilePath,"/",ListOfProtNames[i],"_",ListOfProtNames[j],"_",pot,"_",n,sep=""),
                  row.names = F)
      }
      count <- count + 1
    }
    setTkProgressBar(pb, count+a*total,label=paste( round(count/(2*total)*100+a*50, 0),
                                                "% done"))
  }
}


AllvsAll.EMD <- function(ListOfProtNames,PathToOutput,n,pot)
{
  fileextension <- paste("_",pot,"_",n,sep="")
  
  emdcounter <- 1
  k <- 1
  
  EMD <- numeric(1)
  NameforEMDA <- numeric(1)
  NameforEMDB <- numeric(1)
  
  for(ProtA in ListOfProtNames)
  {
    FLBAA <- read.csv(paste(PathToOutput,ProtA,fileextension,"/",ProtA,"_",ProtA,fileextension,sep=""),header = T)[[1]]
    
    for(ProtB in ListOfProtNames[-c(1:k)])
    {
      FLBAB <- read.csv(paste(PathToOutput,ProtA,fileextension,"/",ProtA,"_",ProtB,fileextension,sep=""),header = T)[[1]]
      
      P <- t(as.matrix(hist(FLBAA , breaks = seq(from = 0,to = max(FLBAA,FLBAB)+0.05, by = 0.05), plot = F)$counts))
      Q <- t(as.matrix(hist(FLBAB , breaks = seq(from = 0,to = max(FLBAA,FLBAB)+0.05, by = 0.05), plot = F)$counts))
      EMD[emdcounter] <- emd2d(Q,P)
      NameforEMDA[emdcounter] <- ProtA
      NameforEMDB[emdcounter] <- ProtB
      emdcounter <- emdcounter + 1
    }
    k <- k+1
  }
  
  write.csv(file=paste(PathToOutput,"ListofEMD_",pot,"_",n,".csv",sep=""),
            x=data.frame(NameforEMDA,NameforEMDB,EMD),row.names = F)
}



AllvsAll.Cluster <- function(PathToOutput,n)
{  
  mydendrogramplot <- function(clust,xlim=NULL,ylim=NULL, title=NULL)
  {
    
    dendrogram <- as.dendrogram(clust)
    dendro.data <- dendro_data(dendrogram)
    
    p <- ggplot() +
      geom_segment(data = dendro.data$segments,
                   aes_string(x = "x", y = "y", xend = "xend", yend = "yend"))+
      theme_dendro()+
      scale_x_continuous(breaks = seq_along(dendro.data$labels$label), 
                         labels = dendro.data$labels$label) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
      ggtitle(title)
    
    if(is.null(xlim) &is.null(ylim))
    {
      p <- p +  coord_cartesian(xlim = xlim, ylim = ylim)
      
    }
    p
  }
  
  pathtoEMDNeg <- paste(PathToOutput,"/ListOfEMD_negative_",n,".csv",sep="")
  pathtoEMDPos <- paste(PathToOutput,"/ListOfEMD_positive_",n,".csv",sep="")
  
  data <- data.frame( read.csv(file= pathtoEMDPos), neg = read.csv(file= pathtoEMDNeg)[,3])
  
  ProtList <- unique(c(as.character(data[,1]),as.character(data[,2])))
  
  matr.Sum <- matrix(0,nrow = NROW(ProtList),ncol = NROW(ProtList), dimnames = list(ProtList,ProtList))
  matr.Pos <- matrix(0,nrow = NROW(ProtList),ncol = NROW(ProtList), dimnames = list(ProtList,ProtList))
  matr.Neg <- matrix(0,nrow = NROW(ProtList),ncol = NROW(ProtList), dimnames = list(ProtList,ProtList))
  matr.Max <- matrix(0,nrow = NROW(ProtList),ncol = NROW(ProtList), dimnames = list(ProtList,ProtList))
  k <- 1
  for(i in 1:NROW(ProtList))
  {
    for(j in (i+1):NROW(ProtList))
    {
      if(k <= NROW(data)){
        matr.Pos[i,j] <- data[k,3]
        matr.Neg[i,j] <- data[k,4]
        matr.Pos[j,i] <- data[k,3]
        matr.Neg[j,i] <- data[k,4]
        matr.Max[i,j] <- max(data[k,3],data[k,4])
        matr.Max[j,i] <- matr.Max[i,j]
        matr.Sum[i,j] <- 1/2*(data[k,3]+data[k,4])
        matr.Sum[j,i] <- matr.Max[i,j]
      }
      k<- k+1
    }
  }
  
  agnes.average.Max <- agnes(x = matr.Max, diss = T,method = "average",keep.diss = F,keep.data = F)
  agnes.average.Neg <- agnes(x = matr.Neg, diss = T,method = "average",keep.diss = F,keep.data = F)
  agnes.average.Pos <- agnes(x = matr.Pos, diss = T,method = "average",keep.diss = F,keep.data = F)
  agnes.average.Sum <- agnes(x = matr.Sum, diss = T,method = "average",keep.diss = F,keep.data = F)
  
  
  mydendrogramplot(agnes.average.Max,title = "UPGMA, Maximum(Neg,Pos)")
  ggsave(filename = paste(PathToOutput,"/Dendrogram_UPGMA_Max.pdf",sep=""),height=7, width = 14)
  
  mydendrogramplot(agnes.average.Neg,title = "UPGMA, only negative")
  ggsave(filename = paste(PathToOutput,"/Dendrogram_UPGMA_Neg.pdf",sep=""),height=7, width = 14)
  
  mydendrogramplot(agnes.average.Pos,title = "UPGMA, only positive")
  ggsave(filename = paste(PathToOutput,"/Dendrogram_UPGMA_Pos.pdf",sep=""),height=7, width = 14)
  
  mydendrogramplot(agnes.average.Sum,title = "UPGMA, Mean")
  ggsave(filename = paste(PathToOutput,"/Dendrogram_UPGMA_Mean.pdf",sep=""),height=7, width = 14)
}

