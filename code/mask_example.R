# this file contains an example of how to run through the image
# processing code. It runs through all 96 images, saves the data frame,
# then builds up an image containing zoomed in, aligned images of all
# identified nuclei in all channels. 

library('EBImage')
source('create_mask.R')

# run through and read in all the images
data <- c()
for (i in 1:96)
{
  id <- sprintf("%05d", i)
  #cat(paste('Processing image', id,'\n'))
  bright <- readImage(paste0('~/Desktop/18511_18855/20170916041443_/AllImages/BRIGHT/',id,'.TIFF'))
  dapi   <- readImage(paste0('~/Desktop/18511_18855/20170916041443_/AllImages/DAPI/',  id,'.TIFF'))
  gfp    <- readImage(paste0('~/Desktop/18511_18855/20170916041443_/AllImages/GFP/',   id,'.TIFF'))
  rfp    <- readImage(paste0('~/Desktop/18511_18855/20170916041443_/AllImages/RFP/',   id,'.TIFF'))

  tmp <- create_mask(bright, dapi, gfp, rfp, id, display=FALSE, size_cutoff=350)
  data <- rbind(data, tmp)
}

write.csv(data, file='imagedata.csv')

# grab all the rows with exactly one nucleus
singles <- data[which(data$total.nuclei==1),]

# make a nice figure. This code is a terrible mess. It creates red,
# green, and blue images for every row/column. I show bg subtracted
# data for red and green channels. 

red.1   <- list()
green.1 <- list()
blue.1  <- list()

red.2   <- list()
green.2 <- list()
blue.2  <- list()

red.3   <- list()
green.3 <- list()
blue.3  <- list()

red.4   <- list()
green.4 <- list()
blue.4  <- list()

for (i in 1:nrow(singles))
{
  bright <- readImage(paste0(singles$id[i],'.1.bright.TIFF'))
  nucs   <- readImage(paste0(singles$id[i],'.1.outline.TIFF'))
  gfps   <- readImage(paste0(singles$id[i],'.1.gfp.TIFF'))
  rfps   <- readImage(paste0(singles$id[i],'.1.rfp.TIFF'))
              
  # convert from images to 3D matrices
  bright <- imageData(rgbImage(bright,bright,bright)) # to make greyscale, red=blue=green
  nucs   <- imageData(nucs)
  gfps   <- imageData(rgbImage(green=gfps))
  rfps   <- imageData(rgbImage(red=rfps))
  
  red.1[[i]]   <- bright[,,1]
  green.1[[i]] <- bright[,,2]
  blue.1[[i]]  <- bright[,,3]
  
  red.2[[i]]   <- nucs[,,1]
  green.2[[i]] <- nucs[,,2]
  blue.2[[i]]  <- nucs[,,3]
  
  red.3[[i]]   <- gfps[,,1]
  green.3[[i]] <- gfps[,,2] - singles$bg.gfp[i] # subtract background
  blue.3[[i]]  <- gfps[,,3]
  
  red.4[[i]]   <- rfps[,,1] - singles$bg.rfp[i] # subtract background
  green.4[[i]] <- rfps[,,2]
  blue.4[[i]]  <- rfps[,,3]
}

reds.1   <- as.Image(do.call('cbind', red.1))
greens.1 <- as.Image(do.call('cbind', green.1))
blues.1  <- as.Image(do.call('cbind', blue.1))

reds.2   <- as.Image(do.call('cbind', red.2))
greens.2 <- as.Image(do.call('cbind', green.2))
blues.2  <- as.Image(do.call('cbind', blue.2))

reds.3   <- as.Image(do.call('cbind', red.3))
greens.3 <- normalize(as.Image(do.call('cbind', green.3))) # normalize to bring scale back up to 0-1 after bg removal
blues.3  <- as.Image(do.call('cbind', blue.3))

reds.4   <- normalize(as.Image(do.call('cbind', red.4))) # normalize to bring scale back up to 0-1 after bg removal
greens.4 <- as.Image(do.call('cbind', green.4))
blues.4  <- as.Image(do.call('cbind', blue.4))

combined <- rgbImage(
  as.Image(rbind(imageData(reds.1),imageData(reds.2),imageData(reds.3),imageData(reds.4))),
  as.Image(rbind(imageData(greens.1),imageData(greens.2),imageData(greens.3),imageData(greens.4))),
  as.Image(rbind(imageData(blues.1),imageData(blues.2),imageData(blues.3),imageData(blues.4))))
  
  
writeImage(combined, file='combined.TIFF')
