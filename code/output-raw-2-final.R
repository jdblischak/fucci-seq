

######## Step 1: sample filtering
## from code/sampleqc.R
## input: data/eset/*.rds
## output: data/eset-raw.rds

## Prepare data
library("Biobase")
fname <- Sys.glob("data/eset/*.rds")
eset <- Reduce(combine, Map(readRDS, fname))
anno <- pData(eset)

## Total mapped reads cutoff
cut_off_reads <- quantile(anno[anno$cell_number == 0,"mapped"], 0.82)
anno$cut_off_reads <- anno$mapped > cut_off_reads

## Unmapped ratio cutoff
anno$unmapped_ratios <- anno$unmapped/anno$umi
cut_off_unmapped <- quantile(anno[anno$cell_number == 0,"unmapped_ratios"], 0.40)
anno$cut_off_unmapped <- anno$unmapped_ratios < cut_off_unmapped

## ERCC percentage cutoff
anno$ercc_percentage <- anno$reads_ercc / anno$mapped
cut_off_ercc <- quantile(anno[anno$cell_number == 0,"ercc_percentage"], 0.20)
anno$cut_off_ercc <- anno$ercc_percentage < cut_off_ercc

## Number of genes detected cutoff
cut_off_genes <- quantile(anno[anno$cell_number == 0,"detect_hs"], 0.80)
anno$cut_off_genes <- anno$detect_hs > cut_off_genes

## Total molecule outlier
## create 3 groups according to cell number
group_3 <- rep("two",dim(anno)[1])
group_3[grep("0", anno$cell_number)] <- "no"
group_3[grep("1", anno$cell_number)] <- "one"
## create data frame
library(dplyr)
data <- anno %>% dplyr::select(experiment:concentration, mapped, molecules)
data <- data.frame(data, group = group_3)

## perform lda to identify outliers
library("MASS")
data_lda <- lda(group ~ concentration + molecules, data = data)
data_lda_p <- predict(data_lda, newdata = data[,c("concentration", "molecules")])$class
data$data_lda_p <- data_lda_p
## identify the outlier
library("tibble")
outliers_lda <- data %>% rownames_to_column("sample_id") %>% filter(cell_number == 1, data_lda_p == "two")
## create filter
anno$molecule_outlier <- row.names(anno) %in% outliers_lda$sample_id


## Read to molecule conversion outlier
## calculate convertion
anno$ercc_conversion <- anno$mol_ercc / anno$reads_ercc
anno$conversion <- anno$mol_hs / anno$reads_hs
## try lda
data$conversion <- anno$conversion
data$ercc_conversion <- anno$ercc_conversion
data_ercc_lda <- lda(group ~ ercc_conversion + conversion, data = data)
data_ercc_lda_p <- predict(data_ercc_lda,  newdata = data[,c("ercc_conversion", "conversion")])$class
data$data_ercc_lda_p <- data_ercc_lda_p
## identify the outlier
outliers_conversion <- data %>% rownames_to_column("sample_id") %>% filter(cell_number == 1, data_ercc_lda_p == "two")
## create filter
anno$conversion_outlier <- row.names(anno) %in% outliers_conversion$sample_id

## Final filter
anno$filter_all <- anno$cell_number == 1 &
  anno$mol_egfp > 0 &
  anno$valid_id &
  anno$cut_off_reads &
  anno$cut_off_unmapped &
  anno$cut_off_ercc &
  anno$cut_off_genes &
  anno$molecule_outlier == "FALSE" &
  anno$conversion_outlier == "FALSE"

## Move filter_all to the last column
anno <- anno[,c(1:42, 44, 43)]

## Update the phenoData slot
## also move filter_all to the last column
varMetadata_updated <- data.frame(labeDescription=
        c(varMetadata(eset)$labelDescription,
         "QC filter: Is the sample an outlier in total molecule count?"),
        stringsAsFactors = F)
varMetadata_updated$labeDescription <- varMetadata_updated$labeDescription[c(1:42, 44,43)]
rownames(varMetadata_updated) <- c(rownames(varMetadata(eset)), "molecule_outlier")[c(1:42,44,43)]

## update the phenoData slot
pData_updated <- new("AnnotatedDataFrame",
                     data=anno, varMetadata=varMetadata_updated)

phenoData(eset) <- pData_updated


## save the updated eset to eset-raw.rds
saveRDS(eset, file = "data/eset-raw.rds")



######## Step 2: gene filtering
## from gene-filter.Rmd
## from 20421 genes to 11093 genes (including ERCC)
## input: data/eset-raw.rds
## output: data/eset-filtered.rds

count_filter <- exprs(eset)

## remove overly expressed genes
which_over_expressed <- which(apply(count_filter, 1, function(x) any(x>(4^6)) ))


## remove lowly expressed genes
cpm <- t(t(count_filter)/pData(eset)$molecules)*(10^6)
which_lowly_expressed <- which(rowMeans(cpm) < 2)
gene_filter <- unique(c(which_lowly_expressed, which_over_expressed))
genes_to_include <- setdiff(1:nrow(count_filter), gene_filter)

## save a filtered data
eset_filtered <- eset[genes_to_include, pData(eset)$filter_all==TRUE]

saveRDS(eset_filtered, file = "data/eset-filtered.rds")



### Step 3: combine intensity data with sequencing data
## from 923 samples in eset-filtered to 888 samples in the final data
## input: data/eset-filtered.rds; data/intensity.rds, 
##        output/images-normalize-anova.Rmd/pdata.adj.rds
## output: data/eset-final.rds

# import data
ints <- readRDS(file="data/intensity.rds")
eset_filtered <- readRDS(file = "data/eset-filtered.rds")
pdata.adj <- readRDS(file = "output/images-normalize-anova.Rmd/pdata.adj.rds")

# identify samples that have data in both eset-filtered.rds and intensity.rds
pdata <- pData(eset_filtered)
pdata$unique <- paste(pdata$image_individual, sprintf("%05d", pdata$image_label), sep="_")
pdata.adj$unique <- paste(pdata.adj$image_individual, sprintf("%05d", pdata.adj$image_label), sep="_")

sample_include_alldata <- intersect(intersect(ints$unique, pdata$unique), pdata.adj$unique)
length(sample_include_alldata)

# combine datasets
ints_combo <- ints[which(ints$unique %in% sample_include_alldata), ]
pdata.adj_combo <- pdata.adj[which(pdata.adj$unique %in% sample_include_alldata),]

eset_combo <- new("ExpressionSet", 
                  exprs = exprs(eset_filtered)[,which(pdata$unique %in% sample_include_alldata)], 
                  phenoData = phenoData(eset_filtered)[which(pdata$unique %in% sample_include_alldata), ], 
                  featureData = featureData(eset_filtered),
		  experimentData = experimentData(eset_filtered))

pdata_combo <- pData(eset_combo)
pdata_combo$unique <- paste(pdata_combo$image_individual, 
                            sprintf("%05d", pdata_combo$image_label), sep="_")

all.equal(ints_combo$unique, pdata_combo$unique)
all.equal(pdata.adj_combo$unique, pdata_combo$unique)


# compute projected cell time
pc.fucci <- prcomp(subset(pdata.adj_combo, 
                          select=c("rfp.median.log10sum.adjust",
                                   "gfp.median.log10sum.adjust")),
                   center = T, scale. = T)
Theta.cart <- pc.fucci$x
library(circular)
Theta.fucci <- coord2rad(Theta.cart)
Theta.fucci <- (2*pi)-as.numeric(Theta.fucci)


# add intensity data to the phenoData slot
pdata_table <- rbind(varMetadata(phenoData(eset_combo)),
                     "mCherry background-corrected intensity (log10sum)",
                     "EGFP background-corrected intensity (log10sum)",
                     "DAPI background-corrected intensity (log10sum)",
                     "mCherry background-corrected intensity adjusted for C1 batch(log10sum)",
                     "EGFP background-corrected intensity adjusted for C1 batch (log10sum)",
                     "DAPI background-corrected intensity adjusted for C1 batch (log10sum)",
                     "nucleus size",
                     "nucleus perimeter",
                     "nucleus eccentricity", 
                     "cell time")
rownames(pdata_table) <- c(rownames(varMetadata(phenoData(eset_combo))),
                           "rfp.median.log10sum",
                           "gfp.median.log10sum",
                           "dapi.median.log10sum",
                           "rfp.median.log10sum.adjust",
                           "gfp.median.log10sum.adjust",
                           "dapi.median.log10sum.adjust",
                           "size", "perimeter", "eccentricity",
                           "theta")


phenoData(eset_combo) <- new("AnnotatedDataFrame", 
       data = data.frame(pData(eset_combo),
                         rfp.median.log10sum=ints_combo$rfp.median.log10sum,
                         gfp.median.log10sum=ints_combo$gfp.median.log10sum,   
                         dapi.median.log10sum=ints_combo$dapi.median.log10sum,
                         rfp.median.log10sum.adjust=pdata.adj_combo$rfp.median.log10sum.adjust,
                         gfp.median.log10sum.adjust=pdata.adj_combo$gfp.median.log10sum.adjust,   
                         dapi.median.log10sum.adjust=pdata.adj_combo$dapi.median.log10sum.adjust,
                         size=ints_combo$size,
                         perimeter=ints_combo$perimeter,
                         eccentricity=ints_combo$eccentricity,
                         theta=Theta.fucci),
       varMetadata = pdata_table)


# save to final eset object
saveRDS(eset_combo, file = "data/eset-final.rds")


