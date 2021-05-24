#' Estimating the intensity ratios across channels within mass range of Gd, the first 155 channel was selected as the reference channel.
#'
#' @param ff The flowframe object of Gd contaminated CyTOF data.
#' @param PercentList The top-n percent of cells selected for estimating Gd ratios, sorted by mean expression of Gd channels.
#'
#' @return The estimated intensity ratios across channels, the 155 channel was selected as the reference channel.
#' @export
#'
#' @examples
#' library(GdClean)
#'
#' PercentList <- seq(5, 100, 5)
#' GdRatios <- estimateGdRatio(ff, percentList)
estimateGdRatio <- function(ff, PercentList) {
  massInfo <- getMassChannels(ff)

  # calculated based on the selected flowframe
  dataTmp <- ff@exprs
  gdData <- dataTmp[, massInfo$massChannelID]
  gdData <- as.data.frame(gdData)
  gdData$mean <- rowMeans(gdData)

  markerName <- ff@parameters@data$name
  gdMarker <- markerName[massInfo$massChannelID]
  PairChannel <- 2:length(massInfo$massChannelID)

  LinearCoeff <- matrix(0, length(PercentList), length(PairChannel))
  for (ii in 1:length(PercentList)) {
    gdData_Sub <- gdData %>% top_frac(n = PercentList[ii] / 100, wt = mean)
    for (jj in 1:length(PairChannel)) {
      data1 <- gdData_Sub[, 1]
      data2 <- gdData_Sub[, PairChannel[jj]]
      dataTmp <- data.frame(baseChannel = data1, selectChannel = data2)

      lmResult <- lm(selectChannel ~ baseChannel, data = dataTmp)
      slope <- lmResult$coefficients[2]
      LinearCoeff[ii, jj] <- slope
    }
  }

  channelName <- paste(gdMarker[2:length(gdMarker)], gdMarker[1], sep = " to ")
  percentName <- paste0(PercentList, "%")

  rownames(LinearCoeff) <- percentName
  colnames(LinearCoeff) <- channelName

  LinearCoeff <- cbind(tmp = 1, LinearCoeff)
  colnames(LinearCoeff)[1] <- paste(gdMarker[1], gdMarker[1], sep = " to ")

  return(LinearCoeff)
}
