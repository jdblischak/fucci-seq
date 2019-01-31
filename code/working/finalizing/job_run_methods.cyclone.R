#library(scran)

# adapted from cyclone github repository
find.markers <- function(id1,id2,id3, training.data,
                         genes.list, thr.frac=0.5) {#function to find markers: compare id1 vs (id2,id3)

  data1<-training.data[which(rownames(training.data) %in% genes.list),
                       which(colnames(training.data) %in% id1)]#
  data2<-training.data[which(rownames(training.data) %in% genes.list),
                       which(colnames(training.data) %in% id2)]#
  data3<-training.data[which(rownames(training.data) %in% genes.list),
                       which(colnames(training.data) %in% id3)]#

  Ngenes<-length(genes.list)
  Nthr1<-ceiling(length(id1)*thr.frac)
  Nthr2<-ceiling(length(id2)*thr.frac)
  Nthr3<-ceiling(length(id3)*thr.frac)

  couple.markers<-matrix(-1,1,3)
  for(i in 1:(Ngenes-1) ){
    for(j in i:Ngenes){
      diff1<-(data1[i,]-data1[j,])
      diff2<-(data2[i,]-data2[j,])
      diff3<-(data3[i,]-data3[j,])


      if(length(diff1[diff1>0])>=Nthr1 & length(diff2[diff2< 0])>=Nthr2 &
         length(diff3[diff3< 0])>=Nthr3){
        couple.markers<-rbind(couple.markers, c(i,j,1))
      }else if(length(diff1[diff1< 0])>=Nthr1 & length(diff2[diff2> 0])>=Nthr2 &
               length(diff3[diff3> 0])>=Nthr3){
        couple.markers<-rbind(couple.markers, c(j,i,1))
      }
    }
  }
  couple.markers<-couple.markers[-1,]
  return(data.frame(Gene.1=genes.list[couple.markers[,1]],
                    Gene.2=genes.list[couple.markers[,2]],
                    Sign=couple.markers[,3]))
}



# use Leng data to find stage markers
df <- readRDS("/project2/gilad/joycehsiao/fucci-seq/data/rnaseq-previous-studies/leng/HumanLengESC.rds")

library(Biobase)
log2cpm <- log2(10^6*t(t(exprs(df)+1)/colSums(exprs(df))))
log2cpm_quant <- t(scale(t(log2cpm)))
pdata <- readRDS("/project2/gilad/joycehsiao/fucci-seq/data/rnaseq-previous-studies/leng/pdata_filtered.rds")
id.G1 <- as.character(pdata$sample_id[pdata$cell_state == "G1"])
id.G2 <- as.character(pdata$sample_id[pdata$cell_state == "G2"])
id.S <- as.character(pdata$sample_id[pdata$cell_state == "S"])
genes.training <- rownames(log2cpm_quant)
training.data <- log2cpm_quant


#Find marker pairs, i.e., pairs of genes that change their relative ranking in phase 1 compared to phase 2 and 3.
#id1 = names of the samples in phase 1
#id2, id3 = names of the samples in phase 2 and 3
#genes.list = list of genes to consider
#thr.frac = fraction of mismatches allowed (default: 0.5)

#load("/project2/gilad/joycehsiao/cyclone/R/pairs_method/core/pairs_functions.RData")
G1.marker.pairs<-find.markers(id1=id.G1,
                              id2=id.S, id3=id.G2,
                              genes.list=genes.training,
                              training.data=training.data, thr.frac=0.5)

S.marker.pairs<-find.markers(id1=id.S,
                             id2=id.G1, id3=id.G2,
                             genes.list=genes.training,thr.frac=0.5)

G2M.marker.pairs<-find.markers(id1=id.G2,
                               id2=id.G1, id3=id.S,
                               genes.list=genes.training,thr.frac=0.5)


pairs.list <- list(G1=G1.marker.pairs[1:2],
                   G2=G2M.marker.pairs[1:2],
                   S=S.marker.pairs[1:2])

saveRDS(pairs.list,
        file="/project2/gilad/joycehsiao/fucci-seq/data/results/finalizing/leng.pairsgene.cyclone.rds")

#assignments <- cyclone(data_withheld$log2cpm.valid, pairs.list)



