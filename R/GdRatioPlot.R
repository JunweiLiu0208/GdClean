#' Comparing the estimated Gd ratios with the of natural Gd abundance ratios
#'
#' @param GdRatios The calculated Gd Ratios from the 'estimateGdRatio'
#''
#' @return The ggplot object of the compared result.
#' @export
#'
#' @examples
#' library(GdClean)
#' p1 <- GdRatioPlot(GdRatios)
#' p1
GdRatioPlot <- function(GdRatios) {
  GdNatureAbundance <- c(0.0020, 0.0218, 0.1480, 0.2047, 0.1565, 0.2484, 0.2186)
  GdNatureTable <- data.frame(Abundance = GdNatureAbundance)
  rownames(GdNatureTable) <- paste0(c(152, 154, 155, 156, 157, 158, 160), "Gd")

  # Plot
  GdRatioName <- colnames(GdRatios)
  plotList <- list()
  for (ii in 2:length(GdRatioName)) {
    idTmp <- GdRatioName[ii]
    idSub <- unlist(str_split(idTmp, "to"))
    GdID_1 <- unlist(regmatches(idSub[1], gregexpr("[[:digit:]]+", idSub[1])))
    GdID_2 <- unlist(regmatches(idSub[2], gregexpr("[[:digit:]]+", idSub[2])))

    GdID_1 <- intersect(GdID_1, c(152, 154, 155, 156, 157, 158, 160))
    GdID_2 <- intersect(GdID_2, c(152, 154, 155, 156, 157, 158, 160))

    GdID_1 <- paste0(GdID_1, "Gd")
    GdID_2 <- paste0(GdID_2, "Gd")

    GdID_1_Abundance <- GdNatureTable[GdID_1, 1]
    GdID_2_Abundance <- GdNatureTable[GdID_2, 1]

    natureRatio <- GdID_1_Abundance / GdID_2_Abundance
    calRatio <- GdRatios[, ii]

    # plot calculated and nature ratios.
    ratioData <- data.frame(calRatio)
    ratioData$percent <- as.numeric(unlist(regmatches(rownames(ratioData), gregexpr("[[:digit:]]+", rownames(ratioData)))))

    p1 <- ggplot(ratioData, aes(x = percent, y = calRatio)) +
      geom_point() +
      stat_smooth(method = gam, formula = y ~ s(x)) +
      geom_hline(
        yintercept = natureRatio, color = "red",
        linetype = "dashed", size = 1
      ) +
      theme_bw() +
      ggtitle(idTmp) +
      labs(x = "Top (n)% of Cells", y = "Coefficients of fitted Linear Models") +
      scale_x_continuous(breaks = c(1:10) * 10) +
      scale_y_continuous(limits = c(0, 3)) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12)
      )
    plotList[[ii - 1]] <- p1
  }

  p2 <- ggarrange(plotlist = plotList, ncol = length(plotList))
  return(p2)
}
