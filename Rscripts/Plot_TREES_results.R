#
#Various plotting function calls to plot TREES simulation results
#

#This path must be set to where you placed the workshop root folder on your computer
setwd(paste("~/Documents/research/projects/TREES-workshop-2020/",sep=""))

source("Rscripts/plotFunctions.R")

#
#Plot leaf growth results
#
subfolder <- "Examples/Brassica/"
fname1 <- "cc_ww_drought_valid_expt"
fname2 <- "cc_lowN"
title <- "Chinese Cabbage Leaf Growth"
plotname <- "Graphics/ccLeaf.pdf"
plotLeafGrowth(subfolder, fname1, fname2, title, plotname, 122)

#
#Plot midday simulation results
#Note: These scripts assume you have already run Aggregate_to_Daily.R 
#
B73 <- read.csv(paste("Examples/Maize/"," B73 _midday",".csv",sep=""),header=TRUE)
title <- "Maize, LIRF 2013"
plotname <- "Graphics/B73.pdf"
plotDaily(B73, title, plotname,151)

pinon <- read.csv(paste("Examples/Pinon/"," ap2012 _midday",".csv",sep=""),header=TRUE)
title <- "SUMO Pinon 2012 - Groundwater Access"
plotname <- "Graphics/pinon.pdf"
plotDaily(pinon, title, plotname, 0)

pinon_noGW <- read.csv(paste("Examples/Pinon/"," ap2012_noGW _midday",".csv",sep=""),header=TRUE)
title <- "SUMO Pinon 2012 - No Groundwater Access"
plotname <- "Graphics/pinon_noGW.pdf"
plotDaily(pinon_noGW, title, plotname, 0)

cc <- read.csv(paste("Examples/Brassica/"," cc_ww_drought_valid_expt _midday",".csv",sep=""),header=TRUE)
title <- "Chinese Cabbage"
plotname <- "Graphics/cc.pdf"
plotDaily(cc, title, plotname,122)

