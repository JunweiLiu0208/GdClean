getGdChannels = function(ff){

  markerPanel = ff@parameters@data
  gdChannelID = str_which(markerPanel$name,'^Gd')
  gdChannelName = markerPanel$name[gdChannelID]
  gdMetalID = as.numeric(unlist(regmatches(gdChannelName, gregexpr("[[:digit:]]+", gdChannelName))))

  gdInfo = data.frame(gdChannelID,gdMetalID)
  rownames(gdInfo) = paste0('Gd',gdInfo$gdMetalID)
  gdInfo = gdInfo %>% arrange(gdMetalID)

  gdChannelNameAll = paste(rownames(gdInfo),collapse = ' ')
  print(paste0(length(gdChannelID),' Gadolinium Channels Selected: ',gdChannelNameAll))

  # different channels
  if(length(gdChannelID) <= 1){
    errorCondition('GdClean can not work with no more than one Gd channel!')
  }else if(length(gdChannelID) == 2){
    print('2 Gd channels detected, the lower expression one will be select as contamination coefficient')
  }else if(length(gdChannelID) >=2){

  }
  return(gdInfo)

}
