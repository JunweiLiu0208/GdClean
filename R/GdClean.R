
GdClean = function(fs,gdRatio,method){

  gdInfo = getGdChannels(fs[[1]])
  gdChannelID = gdInfo$gdChannelID

  for (ii in 1:length(fs)) {

    dataTmp = fs[[ii]]@exprs
    gdDataTmp = dataTmp[,gdChannelID]

    # estimate k
    estimatek = apply(gdDataTmp,1,function(cellData){

      cellK = 0
      switch(method,
             '1DNorm' = {
               obj_1 = function(x){
                 y = sum(abs(cellData - x*gdRatio))
               }
               res1 = optimize(obj_1,interval = c(0,10000))
               res1 = res1$minimum
               cellK = res1
             },
             '2DNorm' = {
               obj_2 = function(x){
                 y = sqrt(sum((cellData - x*gdRatio)^2))
               }
               res2 = optimize(obj_2,interval = c(0,10000))
               res2 = res2$minimum
               cellK = res2
             },
             'min' = {
               k_vector = as.matrix(cellData / gdRatio)
               cellK  = min(k_vector)
             },
             'mean' = {
               k_vector = as.matrix(cellData / gdRatio)
               cellK  = mean(k_vector)
             },
             'median' = {
               k_vector = as.matrix(cellData / gdRatio)
               cellK  = median(k_vector)
             })

      return(cellK)

    })

    # remove Gd contamination
    gdCompensate = gdDataTmp - matrix(estimatek*gdRatio,dim(gdDataTmp)[1],dim(gdDataTmp)[2])
    gdCompensate[gdCompensate < 0] = 0
    gdCompensate = gdCompensate + rnorm(length(gdRatio),1,1)

    # change RawData
    dataTmp[,gdChannelID] = gdCompensate
    fs[[ii]]@exprs = dataTmp

  }

  return(fs)

}
