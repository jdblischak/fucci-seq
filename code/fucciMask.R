#!/usr/bin/env Rscript

#' @title Analyze image from one well
#'
#' @description Generate a mask surrounding nuclei. Take all data channels.
#'            Use DAPI to create a mask, then return a cropped, smaller image containing
#'            each nucleus and writes the cropped, smaller image to files.
#'
#' @param counts gene by sample expression count matrix (G by N).
#' @param libsize_factors numreic vector of the scale factors for library sizes.
#'                        Default to 1 - no adjustment.
#'
#' @param bright
#' @param dapi
#' @param gfp
#' @param rfp
#' @param id
#'  @param control A list with control arguments, including
#'   \code{size_cutoff} default 300.
#'   \code{display} TRUE/FALSE.
#'   \code{write} TRUE/FALSE.
#'
#' @return data.frame contains summary information
#'    \code{id} cell ID
#'    \code{nucleus} number of nuclei
#'    \code{bg.dapi}
#'    \code{bg.gfp}
#'    \code{bg.rfp}
#'    \code{mean.dapi} mean foreground green intensity
#'    \code{mean.gfp} mean foreground green intensity
#'    \code{mean.rfp} mean foreground red intensity
#'    \code{nuc.volume} size of nucleus
#'
#' @example
#' library(EBImage)
#' id <- "00002"
#' dir_images_data_pl="/project2/gilad/fucci-seq/images_curated/18870_18511/"
#' dir_output="/scratch/midway2/joycehsiao/fucci-seq/images_processed/18870_18511/"
#'
#' bright <- readImage(paste0(dir_images_data_pl, "BRIGHT/", id, ".TIFF"))
#' dapi <- readImage(paste0(dir_images_data_pl, "DAPI/", id, ".TIFF"))
#' gfp <- readImage(paste0(dir_images_data_pl, "GFP/", id, ".TIFF"))
#' rfp <- readImage(paste0(dir_images_data_pl, "RFP/", id, ".TIFF"))
#' tmp <- create_mask(bright, dapi, gfp, rfp, id, dir_output = dir_output)

#' @export
create_mask <- function(bright, dapi, gfp, rfp, id,
                        dir_output=NULL,
                        control=list(medianFilterRadius=10,
                                     size_cutoff=200, display=FALSE,
                                     printWholeDAPI=FALSE,
                                     printProcessedImages=FALSE,
                                     displayProcessedImages=FALSE))
{
  library("EBImage")
  library("methods")

  # 1. Work with a normalized image for edge detection:
  # EBImage::normalize() performs linear interpolation of the intensity values to
  # default range (0,1)
  dims <- dim(dapi)
  mask <- normalize(dapi)
  #display(mask, title='Normalized image')

  # 2. Smooth image using 1 round median filtering to improve edge detection:
  # EBImage::medianFilter process an image using Perreault's modern constant-time median
  # filtering algorithm. The main idea to to run through every element in the intensity
  # matrix and replace every element with the median of neighborhood entires. This algorithm
  # defines the neighborhood as a 2*(size+1) by 2*(size+1) square.
  mask <- medianFilter(mask, control$medianFilterRadius)
  #display(mask, title='Median filtered image')

  # 3. Adaptive thresholding
  # EBImage::thresh(): convert to a binary image using a thresholding function.
  # Here we use a 10 pixel width. These could also be user-set parameters
  mask <- thresh(mask, 10, 10, 0.01)
  #display(mask, title='Binary mask')

  # 4. Post-processing filling holes in the image
  # EBImage::fillHull()
  mask <- fillHull(mask)
  #display(mask, title='After fillHull')

  # 5. Post processing smoothing the mask using a erosion/dilation. This removes
  # small features at the border of the mask. Brush width must be odd.
  # EBI::open()
  mask <- opening(mask, makeBrush(9, shape='disc'))
  #display(mask, title='After erosion/dilation')

  # 6. Post processing filling holes the second time. In some images the opening function generated
  # holes, so this will fill them.
  # EBImage::fillHull().
  mask <- fillHull(mask)
  #display(mask, title='After fillHull')

  # 7. Postprocessing labeling nulceus.
  # Label each nucleus with a different integer greyscale value
  # EBImage::bwlabel()
  labeled <- bwlabel(mask)
  nnuc <- max(labeled)

  if (nnuc==0)
  {
    cat(paste("Did not find any nuclei in image", id, "\n"))
    output <- list(id = id,
                   total.nuclei = NA,
                   nuc.volume  = NULL,
                   imageOutput = NULL)
    return(output)

    } else {
  # Remove identified nuclei that are smaller than a user-set size
  # We should think hard about whether this is a necessary step,
  # and it is probably good to visually inspect channels that print
  # this warning
    for (i in 1:nnuc)
    {
      pixels <- which(labeled==i)
      n <- length(pixels)
      if (n < control$size_cutoff)
      {
        cat(paste("Removing an object with size", n, "pixels from image", id, "\n"))
        mask[pixels] <- 0
      }
    }
  }
  labeled <- bwlabel(mask)
  #display(mask, title='After size cutoff')
  nnuc <- max(labeled)

  # if no nuclei was found, return null results
  if (nnuc == 0) {
    output <- list(id = id,
                 total.nuclei = nnuc,
                 nuc.volume  = NULL,
                 imageOutput = NULL)
  }

  # if one or more nuclei were found
  # 8. Grab a 100 pixel image surrounding the nucleus
  # we could consider saving the positions of features as well, so
  # we know if things are in the capture site or not, or how
  # separated features are

  if (nnuc>=1)
  {
    cat(paste("Found", nnuc, "nuclei in image", id, "\n"))

    features <- computeFeatures.moment(labeled)
    start.x <- round(features[,1])-50
    start.y <- round(features[,2])-50
    end.x   <- start.x + 99
    end.y   <- start.y + 99

    # shrink the image if this nucleus is at the image border
    start.x[which(start.x < 1)] <- 1
    start.y[which(start.y < 1)] <- 1
    end.x[which(end.x > dims[1])] <- dims[1]
    end.y[which(end.y > dims[2])] <- dims[2]

    # 9. Print images
    # Show the whole dapi image with identified nuclei highlighted
    if (control$printWholeDAPI == TRUE)
    {
      # print normalized dapi intensities (before processing) in blue
      cells     <- rgbImage(blue=normalize(dapi))
      # overlay nuclei boundary on the normalized intensity image
      # using boundaries detected (labeled)
      segmented <- paintObjects(labeled, cells, col='#ff00ff')
      writeImage(segmented, file=paste0(dir_output,id,'.outlines.TIFF'))
      writeImage(segmented, file=paste0(dir_output,id,'.outlines.png'))
    }

    if (control$printProcessedImages == TRUE) {
      for (i in 1:nnuc)
      {
        name <- paste0(id,'.',i)
        # a zoomed in image of each channel
        this.mask   <- mask[start.x[i]:end.x[i], start.y[i]:end.y[i]]
        this.dapi   <- dapi[start.x[i]:end.x[i], start.y[i]:end.y[i]]
        this.bright <- bright[start.x[i]:end.x[i], start.y[i]:end.y[i]]
        this.rfp    <- rfp[start.x[i]:end.x[i], start.y[i]:end.y[i]]
        this.gfp    <- gfp[start.x[i]:end.x[i], start.y[i]:end.y[i]]


        # zoomed in outlines
        this.label  <- bwlabel(this.mask)
        this.cell   <- rgbImage(blue=normalize(this.dapi))
        segmented   <- paintObjects(this.label, this.cell, col='#ff00ff')

        # write the zoomed images and outlines
        writeImage(segmented,   file=paste0(dir_output, name,'.outline.TIFF'))
        writeImage(this.mask,   file=paste0(dir_output, name,'.mask.TIFF'))
        writeImage(this.bright, file=paste0(dir_output, name,'.bright.TIFF'))
        writeImage(this.dapi,   file=paste0(dir_output, name,'.dapi.TIFF'))
        writeImage(this.rfp,    file=paste0(dir_output, name,'.rfp.TIFF'))
        writeImage(this.gfp,    file=paste0(dir_output, name,'.gfp.TIFF'))

        # show them if you want
        if (control$displayProcessedImages == TRUE)
        {
          display(segmented)
          display(this.bright)
          display(this.dapi)
          display(this.rfp)
          display(this.gfp)
        }
      }
    }

    # The area of the nucleus, in pixels. Other metrics are available
    size <- computeFeatures.shape(labeled)[,1]

    output <- list(id = id,
                 total.nuclei = nnuc,
                 nuc.volume  = size,
                 imageOutput = list(label = labeled,
                                    dapi = dapi,
                                    rfp = rfp,
                                    gfp = gfp))
   }

  # grab the mean fluorescence within identified nuclei
  # other metrics are available if we want them
  #  mean.dapi <- computeFeatures.basic(labeled, dapi)[,1]
  #  mean.rfp  <- computeFeatures.basic(labeled, rfp)[,1]
  #  mean.gfp  <- computeFeatures.basic(labeled, gfp)[,1]

  return( output )
}
