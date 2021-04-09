
estimateGdRatio = function(ff,PercentList){

 gdInfo = getGdChannels(ff)

# calculated based on the selected flowframe
 dataTmp = ff@exprs
 gdData = dataTmp[,gdInfo$gdChannelID]
 gdData = as.data.frame(gdData)
 gdData$mean = rowMeans(gdData)

 markerName = ff@parameters@data$name
 gdMarker = markerName[gdInfo$gdChannelID]

 PairChannel = 2:length(gdInfo$gdChannelID)

 LinearCoeff = matrix(0,length(PercentList),length(PairChannel))
 for (ii in 1:length(PercentList)) {

    gdData_Sub = gdData %>% top_frac(n = PercentList[ii] / 100,wt = mean)
    for (jj in 1:length(PairChannel)) {

      data1 = gdData_Sub[,1]
      data2 = gdData_Sub[,PairChannel[jj]]
      dataTmp = data.frame(baseChannel = data1,selectChannel = data2)

      lmResult = lm(selectChannel ~ baseChannel,data = dataTmp)
      slope = lmResult$coefficients[2]
      LinearCoeff[ii,jj] = slope

    }
 }

channelName = paste(gdMarker[2:length(gdMarker)],gdMarker[1],sep = ' to ')
percentName = paste0(PercentList,'%')

rownames(LinearCoeff) = percentName
colnames(LinearCoeff) = channelName

LinearCoeff = cbind(tmp = 1,LinearCoeff)
colnames(LinearCoeff)[1] = paste(gdMarker[1],gdMarker[1],sep = ' to ')

return(LinearCoeff)

}
