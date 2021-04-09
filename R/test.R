rm(list = ls())

library(flowCore)
library(stringr)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(mgcv)

# load fcs files
fcs.dir = '../RawData/'
frames = lapply(dir(fcs.dir,full.names = T),read.FCS,transformation = F,alter.names = T)
names(frames) <- basename(dir(fcs.dir))
fs = as(frames,'flowSet')

# Plot Gd Correlation Heatmap
p1 = GdCorrHeatmap(fs)
p1

# Estimate Gd Ratios
PercentList = seq(5,100,5)
ff = fs[[1]]
GdRatios = estimateGdRatio(ff,percentList)

# plot Gd Ratios
p2 = GdRatioPlot(GdRatios)
p2

# Gd Clean
gdRatio = GdRatios['5%',]
method = '1DNorm'
fs_Clean = GdClean(fs,gdRatio = gdRatio,method = method)

#
p3 = GdCorrHeatmap(fs_Clean)
p3





