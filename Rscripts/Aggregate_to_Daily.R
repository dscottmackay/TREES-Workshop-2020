#
#TREES output aggregation to daily
#  This script is used first to convert half-hourly simulation output to daily
#

####################################
## set your working directory and load functions
####################################
setwd(paste("~/Documents/research/projects/TREES-workshop-2020/",sep=""))
source("Rscripts/TREES_aggregation_function.R")

#
#run these lines to create daily files for the four half-hourly simulations
#
#use these lines for chinese cabbage
subfolder <- "Examples/Brassica/"
fname <- "cc_ww_drought_valid_expt"
driver<- read.table(paste(subfolder,"cc_ww_drivers_drought_validation",".txt",sep=""), header=TRUE)
Ksat <- 4.0/(-0.25+0.9) #whole-plant saturated hydraulic conductance
simulation<-read.table(paste(subfolder,fname,".sim",sep=""),header=TRUE)
computeDaily(subfolder, fname, simulation, driver, Ksat)


#use these lines for pinon - bedrock water access
subfolder <- "Examples/Pinon/"
fname <- "ap2012"
driver<- read.table(paste(subfolder,"ap2012-met",".txt",sep=""), header=TRUE)
Ksat <- 1.5/(-0.93+2.13) #whole-plant saturated hydraulic conductance
simulation<-read.table(paste(subfolder,fname,".sim",sep=""),header=TRUE)
computeDaily(subfolder, fname, simulation, driver, Ksat)

#pinon - no bedrock water access
subfolder <- "Examples/Pinon/"
fname <- "ap2012_noGW"
driver<- read.table(paste(subfolder,"ap2012-met",".txt",sep=""), header=TRUE)
Ksat <- 1.5/(-0.93+2.13) #whole-plant saturated hydraulic conductance
simulation<-read.table(paste(subfolder,fname,".sim",sep=""),header=TRUE)
computeDaily(subfolder, fname, simulation, driver, Ksat)


#use these lines for maize
subfolder <- "Examples/Maize/"
fname <- "B73"
driver<- read.table(paste(subfolder,"LIRFcorn_2013_grow_irrigated",".txt",sep=""), header=TRUE)
Ksat <- 4.22/(-0.3+1.5) #whole-plant saturated hydraulic conductance
simulation<-read.table(paste(subfolder,fname,".sim",sep=""),header=TRUE)
computeDaily(subfolder, fname, simulation, driver, Ksat)

#use these lines for maize
subfolder <- "Examples/Maize/"
fname <- "MO18W"
driver<- read.table(paste(subfolder,"LIRFcorn_2013_grow_irrigated",".txt",sep=""), header=TRUE)
Ksat <- 4.22/(-0.3+2.03) #whole-plant saturated hydraulic conductance
simulation<-read.table(paste(subfolder,fname,".sim",sep=""),header=TRUE)
computeDaily(subfolder, fname, simulation, driver, Ksat)

#use these lines for maize
subfolder <- "Examples/Maize/"
fname <- "CML103"
driver<- read.table(paste(subfolder,"LIRFcorn_2013_grow_irrigated",".txt",sep=""), header=TRUE)
Ksat <- 4.22/(-0.3+1.57) #whole-plant saturated hydraulic conductance
simulation<-read.table(paste(subfolder,fname,".sim",sep=""),header=TRUE)
computeDaily(subfolder, fname, simulation, driver, Ksat)


