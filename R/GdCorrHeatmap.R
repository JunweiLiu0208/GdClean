#' Calculate the Pearson correlation coefficients across Gd channels in each CyTOF data in the flowSet object.
#'
#' @param fs The input flowSet object
#'
#' @return The output of pheatmap object
#' @export
#'
#' @examples
#' library(GdClean)
#' p1 <- GdCorrHeatmap(fs)
#' p1
GdCorrHeatmap <- function(fs) {
  massInfo <- getMassChannels(fs[[1]])

  channelPairs <- combn(massInfo$massChannelID, 2)
  corrHeatmap <- matrix(0, length(fs), dim(channelPairs)[2])
  sampleID <- sampleNames(fs)

  for (ii in 1:length(fs)) {
    sampleData <- fs[[ii]]@exprs
    corrTmp <- apply(channelPairs, 2, function(x) {
      y <- cor(sampleData[, x[1]], sampleData[, x[2]])
      y
    })
    corrHeatmap[ii, ] <- corrTmp
  }

  markerName <- fs[[1]]@parameters@data$name
  name1 <- markerName[channelPairs[1, ]]
  name2 <- markerName[channelPairs[2, ]]
  pairName <- paste(name1, name2, sep = "_")

  rownames(corrHeatmap) <- sampleID
  colnames(corrHeatmap) <- pairName

  breaksList <- seq(-1, 1, by = 0.1)
  colorList <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(length(breaksList))

  p1 <- pheatmap(corrHeatmap,
    cluster_rows = F, cluster_cols = F,
    color = colorList,
    breaks = breaksList, angle_col = 90,
    display_numbers = round(corrHeatmap, 2)
  )
  p1

  return(p1)
}
