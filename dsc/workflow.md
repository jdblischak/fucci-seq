

dir_main <- "/project2/gilad/joycehsiao/fucci-seq"

1. import data from `dir_main/data/eset-final.rds`

2. prepare data
   input: `dir_main/data/eset-final.rds`,
          `dir_main/output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds`
   output: `dir_main/data/results` 
   code: makedata_supervised.R
   steps: eset 
