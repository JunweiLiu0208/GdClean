#' Return the information of channels with Gd isotope mass in CyTOF data.
#'
#' @param ff flowFrame object
#'
#' @return dataFrame of channel information within the mass range of Gd isotopes.
#' @export
#'
#' @examples
#' library(GdClean)
#' massInfo <- getMassChannels(ff)
getMassChannels <- function(ff) {
  markerPanel <- ff@parameters@data
  massID <- list("155", "156", "157", "158", "160")
  massChannelID <- sapply(massID, function(x) {
    out <- str_which(markerPanel$name, x)
    out
  })

  massChannelName <- markerPanel$name[massChannelID]
  massMarkerName <- markerPanel$desc[massChannelID]

  massInfo <- data.frame(massChannelID, massChannelName, massMarkerName)
  rownames(massInfo) <- massID

  ChannelNameAll <- paste(massInfo$massMarkerName, collapse = " ")
  print(paste0(length(massChannelID), " Channels Selected: ", ChannelNameAll))

  return(massInfo)
}
