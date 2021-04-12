
#' Removing Gd contamination in CyTOF data.
#'
#' @param fs The flowSet object of contaminated CyTOF files.
#' @param gdRatio The selected intensity ratios across Gd channels, default with intensity ratios of natural abundance
#' @param method The method used for calculating contamination coefficients for each cell, default with '1DNorm', available options are '1DNorm','2DNorm','Min','Mean','Median'
#'
#' @return The flowSet object of cleaned CyTOF files.
#' @export
#'
#' @examples
#' library(GdClean)
#'
#' fcs.dir <- "../fcsData/"
#' frames <- lapply(dir(fcs.dir, full.names = T), read.FCS, transformation = F, alter.names = T)
#' names(frames) <- basename(dir(fcs.dir))
#' fs <- as(frames, "flowSet")
#'
#' PercentList <- seq(5, 100, 5)
#' GdRatios <- estimateGdRatio(ff[[1]], percentList)
#'
#' gdRatio <- GdRatios["5%", ]
#' method <- "1DNorm"
#' fs_Clean <- GdClean(fs, gdRatio = gdRatio, method = method)
GdClean <- function(fs, gdRatio = NULL, method = "1DNorm") {

  gdInfo <- getGdChannels(fs[[1]])
  gdChannelID <- gdInfo$gdChannelID

  GdNatureAbundance <- c(0.0020, 0.0218, 0.1480, 0.2047, 0.1565, 0.2484, 0.2186)
  GdNatureTable <- data.frame(Abundance = GdNatureAbundance)
  rownames(GdNatureTable) <- paste0(c(152, 154, 155, 156, 157, 158, 160), "Gd")

  naturalRatio <- GdNatureTable[rownames(gdInfo), "Abundance"]
  naturalRatio <- naturalRatio / naturalRatio[1]

  if (is.null(gdRatio)) {
    gdRatio <- naturalRatio
  }

  print(paste0("Gd Ratios = ", paste(gdRatio, collapse = ", ")))

  if (length(gdChannelID) == 2) {
    for (ii in 1:length(fs)) {
      dataTmp <- fs[[ii]]@exprs
      gdDataTmp <- dataTmp[, gdChannelID]

      # estimate k
      # only two Gd channels select the lower expression as contamination coefficient
      estimatek <- apply(gdDataTmp, 1, min)

      # remove Gd contamination
      gdRatio = matrix(gdRatio,nrow = 1)
      gdNoiseData = gdRatio[rep(1,length(estimatek)), ] * estimatek

      gdCompensate <- gdDataTmp - gdNoiseData
      gdCompensate[gdCompensate < 0] <- 0 # remove negative values

      noiseData = matrix(rnorm(dim(gdDataTmp)[1]*dim(gdDataTmp)[2], 1, 1),
                         dim(gdDataTmp)[1], dim(gdDataTmp)[2])
      gdCompensate <- gdCompensate + noiseData # add noise

      # replace expression data
      dataTmp[, gdChannelID] <- gdCompensate
      fs[[ii]]@exprs <- dataTmp
    }
  } else {
    for (ii in 1:length(fs)) {
      dataTmp <- fs[[ii]]@exprs
      gdDataTmp <- dataTmp[, gdChannelID]

      # estimate k
      estimatek <- apply(gdDataTmp, 1, function(cellData) {
        cellK <- 0
        switch(method,
          "1DNorm" = {
            obj_1 <- function(x) {
              y <- sum(abs(cellData - x * gdRatio))
            }
            res1 <- optimize(obj_1, interval = c(0, 10000))
            res1 <- res1$minimum
            cellK <- res1
          },
          "2DNorm" = {
            obj_2 <- function(x) {
              y <- sqrt(sum((cellData - x * gdRatio)^2))
            }
            res2 <- optimize(obj_2, interval = c(0, 10000))
            res2 <- res2$minimum
            cellK <- res2
          },
          "Min" = {
            k_vector <- as.matrix(cellData / gdRatio)
            cellK <- min(k_vector)
          },
          "Mean" = {
            k_vector <- as.matrix(cellData / gdRatio)
            cellK <- mean(k_vector)
          },
          "Median" = {
            k_vector <- as.matrix(cellData / gdRatio)
            cellK <- median(k_vector)
          }
        )
        return(cellK)
      })

      ##
      gdRatio = matrix(gdRatio,nrow = 1)
      gdNoiseData = gdRatio[rep(1,length(estimatek)), ] * estimatek

      gdCompensate <- gdDataTmp - gdNoiseData
      gdCompensate[gdCompensate < 0] <- 0 # remove negative values

      noiseData = matrix(rnorm(dim(gdDataTmp)[1]*dim(gdDataTmp)[2], 1, 1),
                         dim(gdDataTmp)[1], dim(gdDataTmp)[2])
      gdCompensate <- gdCompensate + noiseData # add noise

      # replace expression data
      dataTmp[, gdChannelID] <- gdCompensate
      fs[[ii]]@exprs <- dataTmp
    }
  }
  return(fs)
}
