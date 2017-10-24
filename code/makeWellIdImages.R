#!/usr/bin/env Rscript

# Description:
#   This script was used to create images of well number ID 
# from 00001 to 00096. The script may be later adpated to print
# well ID in the format of column and row numbers.
 
dir_wellIdImages <- "/project2/gilad/fucci-seq/images_wellid/"

for (i in 1:96)
{
  id <- sprintf("%05d", i)

  tiff(paste0(dir_wellIdImages, id, ".tiff"), width = 100, height = 100)
  par(mfrow=c(1,1), mar=rep(0,4), oma=rep(0,4))
  plot(c(1:4),c(1:4), pch = "", axes = FALSE, ann = FALSE)
  text(x=2.5, y = 2.5, labels=id, adj = c(.5, .5), cex = 1.5, font = 2)
  dev.off()
}
