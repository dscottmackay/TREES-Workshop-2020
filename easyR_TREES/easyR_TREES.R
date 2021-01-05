###############################################################################################################
###############################################################################################################
## Run examples:
## easyR_TREES 
## version 1.1
## 28122020
## Mitchell Hitchcock
## mchitchc@buffalo.edu
## 
## Description:
## An example script for modifying parameters, running TREES,  making common graphs, and other
## useful functions to make running multiple TREES simulations easier.
##
###############################################################################################################
###############################################################################################################


###############################################################################################################
###############################################################################################################
############################### Examples of Running and plotting easyR_TREES ##################################
###############################################################################################################
###############################################################################################################




##set working directory to where you downloaded off of github. Example if you downloaded it to your desktop:
#setwd("~/Desktop/TREES-workshop-2020")


#Run these two lines after setting your working directory.
source("easyR_TREES/Functions_easyR_TREES.R")
source("easyR_TREES/Graph_easyR_TREES.R")




###############################################################################################################
########################################## Brassica: Chinese Cabbage ##########################################
###############################################################################################################

#Run the TREES model with specified parameters. See argument comments for explanations.
Run_easyR_TREES(Result_Dest="easyR_TREES/Outputs",#Your destination for TREES outputs.
                Init_param="Examples/Brassica/cc.p",#Name and location of initial parameter file.
                Which_Parameter="microbiomeScalar",#Parameter to modify (MUST MATCH PARAMETER NAME IN .p file exactly).
                New_values=c(100,200,250,300,500), #A vector of the new values for the specified parameter.
                Use_Gamma=FALSE,#T or F use the growth gamma function.
                Itter=NULL,#Number of iterations for running over the gamma.
                Drivers="cc_ww_drivers_drought_validation",#Vector of names of drivers to use.
                Driver_loc="Examples/Brassica",#Location of drivers if not in working directory.
                Init_Infile="Examples/Brassica/in_cc_ww.txt",#Name and location of initial infile to be modified.
                param_mod = "Examples/Brassica/param_mod",#Name and location of param_mod file.
                covfile="Examples/Brassica/covfile1",#Name and location of covfile.
                N_cores=6,#Number of cores to use for parallel runs. WARNING! This should be at least 1 less than total cores
                TREES_loc="Model-Code/trees3",#Location of *compiled* TREES model.
                show_TREES_cout=TRUE#FALSE will hide the trees console readout.
)


#Generate plots of the leaf areas
Leaf_Area_Plot(Result_Dest="easyR_TREES/Outputs",#Location you set for TREES outputs.
               Which_Parameter="microbiomeScalar",#Parameter that was modified.
               New_values=c(100,200,250,300,500),#Values parameter was changed to.
               Use_Gamma=FALSE,#Set to TRUE if Gamma was used.
               Itter=NULL,#Number of ittereations if Gamma was used.
               Drivers="cc_ww_drivers_drought_validation",#Vector of the drivers used.
               Figure_title=NULL,#Optional title for graphs.
               Compare_to=NULL,#Optional data.frame of experimental data to compare to.
               Smooth_on=FALSE#False will create a boxplot at each day and smooth will create a smooth line with geom_smooth 
                              ##(only applies to exp data or if gamma was on)
)


#Generate plots of 
Sim_Plot(Result_Dest="easyR_TREES/Outputs",
         Which_Sim = "R_total",#Name of the column in the .sim file that you wish to plot (must match exactly).
         Which_Parameter="microbiomeScalar",
         New_values=c(100,200,250,300,500),
         Use_Gamma=FALSE,
         Itter=NULL,
         Drivers="cc_ww_drivers_drought_validation",
         Figure_title=NULL,
         Compare_to=NULL,
         Smooth_on=TRUE
         
)



###############################################################################################################
######################################### Wisconsin Fast Plant Example ########################################
###############################################################################################################

#Run over each driver which corresponds to experimental blocks or different locations.

Run_easyR_TREES(Result_Dest="easyR_TREES/Outputs",
                Init_param="Examples/FP/paramV10.p",
                Which_Parameter="microbiomeScalar",
                New_values=c(10,15,20), 
                Use_Gamma=FALSE,
                Itter=NULL,
                Drivers=c("FP_Drivers_BLK_1",
                          "FP_Drivers_BLK_2",
                          "FP_Drivers_BLK_3",
                          "FP_Drivers_BLK_4",
                          "FP_Drivers_BLK_5",
                          "FP_Drivers_BLK_6"),
                Driver_loc="Examples/FP",
                Init_Infile="Examples/FP/in_FP_by_blk.txt",
                param_mod = "Examples/FP/param_mod",
                covfile="Examples/FP/covfile1",
                N_cores=6,#Remember to change this for your computer
                TREES_loc="Model-Code/trees3",
                show_TREES_cout=TRUE
)


#Check the leaf area (keep in mind FP grows leaves in an odd way so the results are actually projected leaf area)
Leaf_Area_Plot(Result_Dest="easyR_TREES/Outputs",
               Which_Parameter="microbiomeScalar",
               New_values=c(10,15,20),
               Use_Gamma=FALSE,Itter=NULL,
               Drivers=c("FP_Drivers_BLK_1",
                         "FP_Drivers_BLK_2",
                         "FP_Drivers_BLK_3",
                         "FP_Drivers_BLK_4",
                         "FP_Drivers_BLK_5",
                         "FP_Drivers_BLK_6"),
               Figure_title=NULL,
               Compare_to=NULL,
               Smooth_on=FALSE
)




###############################################################################################################
############################# Example of experimental data formatting #########################################
###############################################################################################################
#Experimental data must have four columns named exactly as Jday, Leaf_Area, Leaf_number, and block. 
#Block is ordered in the same way that the drivers were created.
#For example, block 1 is the same as the first driver used in the TREES model.
#
#Below is an example for FP experimental data
#FP was treated as a single leaf equivalent to projected area due to its growth behavior and reported data so 
#Leaf_number in this example is 1.
Exp_data<- data.frame() 
for(i in 1:6){# 1:6 represent the blocks which correspond to the drivers used in the same order.
  
  Exp_data<-rbind(Exp_data,
                  read_excel("easyR_TREES/time_course.xlsx", sheet = "data") %>% #raw experimental data
                    dplyr::select(sampling_day_actual,sampling_time_actual,block,soil_trt,contains("avg")) %>% #Select relevant columns
                    filter(soil_trt=="ATM_BLANK",block==i) %>% #"ATM_BLANK" is no microbial treatment Or "SBC_OLD" for microbial treatment
                    dplyr::select(contains("avg_PA_day")) %>%
                    rename_all(funs(c(1,3,4,6,8,9,10,11,12,13,14,16)))%>% #Rename the columns for each day to number (because data is wide)
                    pivot_longer(cols=colnames(.),names_to = "Day",values_to = "Leaf_Area") %>% #Convert the wide data to long
                    mutate(Day=as.double(.$Day)) %>% 
                    arrange(Day) %>% 
                    mutate(Day=Day+332) %>% #Day 1 is Jday 333 so adding 332 to all the days gets Jday
                    rename(Jday=Day) %>% 
                    mutate(Leaf_number=1,block=i) %>% 
                    filter(!is.na(Leaf_Area))
  )
}  


head(Exp_data)#This is what your experimental data should look like 
##When comparing other variables in Sim_Plot() it will be similar but replacing Leaf_Area 
#with the variable name and eliminating Leaf_number.


#Example plotting with experimental data
Leaf_Area_Plot(Result_Dest="easyR_TREES/Outputs",
               Which_Parameter="microbiomeScalar",
               New_values=c(10,15,20),
               Use_Gamma=FALSE,Itter=NULL,
               Drivers=c("FP_Drivers_BLK_1",
                         "FP_Drivers_BLK_2",
                         "FP_Drivers_BLK_3",
                         "FP_Drivers_BLK_4",
                         "FP_Drivers_BLK_5",
                         "FP_Drivers_BLK_6"),
               Figure_title=NULL,
               Compare_to=Exp_data,#This is where the experimental data data.frame is input
               Smooth_on=FALSE
)






###############################################################################################################
##################################### Example of using Gamma ##################################################
###############################################################################################################
#This is an example using the stochastic leaf growth in the leaf growth module. Since the model will sample
# from a gamma distribution to get an idea of the results distribution you need to run the model through many
# iterations. Below is an example of how to do this using the easyR_TREES script. This is unrealistic as it 
# would likely take thousands of iterations to get an accurate distribution but this version is safe for most
# regular computers to handle. This will run 2 different parameters over 6 drivers for 5 itterations. this will
# make for a total of 60  simulations. 
Run_easyR_TREES(Result_Dest="easyR_TREES/Outputs",
                Init_param="Examples/FP/paramV10.p",
                Which_Parameter="microbiomeScalar",
                New_values=c(15), #Each new value will run the number of iterations over each block.
                Use_Gamma=TRUE,#T or F use the growth gamma function
                Itter=c(1:5),#number of iterations for running over the gamma
                Drivers=c("FP_Drivers_BLK_1",
                          "FP_Drivers_BLK_2",
                          "FP_Drivers_BLK_3",
                          "FP_Drivers_BLK_4",
                          "FP_Drivers_BLK_5",
                          "FP_Drivers_BLK_6"),
                Driver_loc="Examples/FP",
                Init_Infile="Examples/FP/in_FP_by_blk.txt",
                param_mod = "Examples/FP/param_mod",
                covfile="Examples/FP/covfile1",
                N_cores=6,#Remember to change this for your computer
                TREES_loc="Model-Code/trees3",
                show_TREES_cout=FALSE
)



Leaf_Area_Plot(Result_Dest="easyR_TREES/Outputs",
               Which_Parameter="microbiomeScalar",
               New_values=c(15),
               Use_Gamma=TRUE,Itter=c(1:5),
               Drivers=c("FP_Drivers_BLK_1",
                         "FP_Drivers_BLK_2",
                         "FP_Drivers_BLK_3",
                         "FP_Drivers_BLK_4",
                         "FP_Drivers_BLK_5",
                         "FP_Drivers_BLK_6"),
               Figure_title="Leaf area dist for low SMB scalar",
               Compare_to=Exp_data,
               Smooth_on=FALSE
)



Sim_Plot(Result_Dest="easyR_TREES/Outputs",
         Which_Sim = "PlantNstat",
         Which_Parameter="microbiomeScalar",
         New_values=c(15),
         Use_Gamma=TRUE,Itter=c(1:5),
         Drivers=c("FP_Drivers_BLK_1",
                   "FP_Drivers_BLK_2",
                   "FP_Drivers_BLK_3",
                   "FP_Drivers_BLK_4",
                   "FP_Drivers_BLK_5",
                   "FP_Drivers_BLK_6"),
         Figure_title=NULL,
         Compare_to=NULL,
         Smooth_on=FALSE
         
)




###############################################################################################################
###################################### Try your own data! #####################################################
###############################################################################################################
# If your own data isn't available you could practice with some of the examples not used here

Run_easyR_TREES(Result_Dest="easyR_TREES/Outputs",#Your destination for TREES outputs.
                Init_param="",#Name and location of initial parameter file.
                Which_Parameter="",#Parameter to modify (MUST MATCH PARAMETER NAME IN .p file exactly).
                New_values=c(), #A vector of the new values for the specified parameter.
                Use_Gamma=FALSE,#T or F use the growth gamma function.
                Itter=NULL,#Number of iterations for running over the gamma.
                Drivers="",#Vector of names of drivers to use.
                Driver_loc="",#Location of drivers if not in working directory.
                Init_Infile="",#Name and location of initial infile to be modified.
                param_mod = "",#Name and location of param_mod file.
                covfile="",#Name and location of covfile.
                N_cores=1,#Number of cores to use for parallel runs. WARNING! This should be at least 1 less than total cores
                TREES_loc="Model-Code/trees3",#Location of *compiled* TREES model.
                show_TREES_cout=TRUE#FALSE will hide the trees console readout.
)









###############################################################################################################
###############################################################################################################
###############################################################################################################



