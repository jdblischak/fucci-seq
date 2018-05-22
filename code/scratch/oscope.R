---
  
  ## Oscope
  
library("Oscope")

Sizes <- MedianNorm(counts)
DataNorm <- GetNormalizedMat(counts, Sizes)
MV <- CalcMV(Data = counts, Sizes = Sizes, Plot = TRUE, MeanCutLow = 1)
DataSubset <- DataNorm[MV$GeneToUse,]

MV2 <- CalcMV(Data = DataSubset, Sizes = NULL, NormData = TRUE, MeanCutLow = 1)
DataSubset2 <- DataNorm[MV2$GeneToUse,]

DataInput <- NormForSine(DataSubset)
SineRes <- OscopeSine(DataInput, parallel = TRUE)
KMRes <- OscopeKM(SineRes, maxK = 10)
ToRM <- FlagCluster(SineRes, KMRes, DataInput)
KMResUse <- KMRes[-ToRM$FlagID]
ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput, NCThre = 100, parallel=TRUE)



