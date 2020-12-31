###############################################################################################################
###############################################################################################################
## Main functions for:
## easyR_TREES 
## version 1.1
## 28122020
## Mitchell Hitchcock
## mchithc@buffalo.edu
## 
## Description:
## A script for modifying parameters, modifying input files, and running TREES.
## Modify at own risk.
###############################################################################################################
###############################################################################################################


##################################Load required packages#######################################################


# List of required packages
pkgli<-c("dplyr","tidyr","readr","readxl","parallel","foreach","doParallel")



## Install missing packages. WARNING!! this can fail one some systems see example snippet following this one
for(pa in pkgli){
  if (!require(pa,character.only=TRUE)){
    install.packages(pa,character.only=TRUE)
  }else{}
  library(pa,character.only=TRUE)
}

rm(pkgli,pa)

###############################################################################################################
##Below is for CCR but also serves as an example of how to modify the above lines as needed.
##This is only ever and issue for packages that are not already installed and if you do not have
## sufficient privileges. 
# for(pa in pkgli){
#   if (!require(pa,character.only=TRUE)){
#     install.packages(pa,character.only=TRUE,repos="https://repo.miserver.it.umich.edu/cran/"$
#   }else{}
#   library(pa,character.only=TRUE,lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
# }
###############################################################################################################


###############################################################################################################
###############################################################################################################
####################################### Modifying Run Functions ###############################################
###############################################################################################################
###############################################################################################################

########################### Function to Modify Trees param.p file #############################################
## The parameter text only has to have a unique portion of the parameter line to be modified
## WARNING!! if your Parameter input has multiple matches ALL of them will be modified
## so it is always better to be more specific with your parameter name.
ChangeParameter<-function(InParamFile="param.p",Parameter="none",NewValue="0"){
  NewValue=paste(NewValue,"")
  ParmFile<-readLines(con=InParamFile)
  ParmFile[grep(Parameter,ParmFile)]<-ParmFile[grep(Parameter,ParmFile)] %>% 
    sub(pattern ='^\\S+ ',replacement=NewValue)
  #return(ParmFile[grep(Parameter,ParmFile)])
  In_param_temp<-tempfile()
  writeLines(ParmFile,In_param_temp)
  return(In_param_temp)
}
###############################################################################################################

########################### Function to Modify the infile for TREES ###########################################
##This is how you will make sure your saved outputs are unique and don't overwrite other runs
##driver is the name of your driver file, OutputPrefix will define both the folder name and the
##unique prefix for each of the runs, and InFileLoc is the location of the original infile
ChangeInFile<-function(driver=NULL,OutputPrefix="outputs/pref",ParamFile="param.p",InFileLoc="in_file.txt",param_mod_loc="param_mod",covfile_loc="covfile1"){
  InFileTREES<-readLines(con=InFileLoc)
  
  if(!is.null(driver)){
    InFileTREES[1]<-driver
  }else{}
  InFileTREES[3]<-ParamFile
  InFileTREES[4]<-param_mod_loc
  InFileTREES[5]<-OutputPrefix
  InFileTREES[7]<-covfile_loc
  
  In_file_temp<-tempfile()
  writeLines(InFileTREES,In_file_temp)
  return(In_file_temp)
}
###############################################################################################################

###############################################################################################################
###############################################################################################################
###############################################################################################################


###############################################################################################################
###############################################################################################################
########################################## Main Run Function ##################################################
###############################################################################################################
###############################################################################################################



Run_easyR_TREES<-function(Result_Dest=getwd(),#Your destination for TREES outputs.
                          Init_param=paste(getwd(),"param.p",sep=""),#Name and location of initial parameter file.
                          Which_Parameter="none",#Parameter to modify (MUST MATCH PARAMETER NAME IN .p file exactly).
                          New_values=NULL, #A vector of the new values for the specified parameter.
                          Use_Gamma=FALSE,#T or F use the growth gamma function.
                          Itter=NULL,#Number of iterations for running over the gamma.
                          Drivers=c("Driver1"),#Vector of names (as character string) of drivers to use.
                          Driver_loc=getwd(),#Location of drivers if not in working directory.
                          Init_Infile="in_FP_by_blk.txt",#Name and location of initial infile to be modified.
                          param_mod="param_mod",#Name and location of param_mod file.
                          covfile="covfile1",#Name and location of covfile.
                          N_cores=2,#Number of cores to use for parrallel runs.
                          TREES_loc="./TREES_model_3.1.2/trees3",#Location of *compiled* TREES model.
                          show_TREES_cout=F#FALSE will hide the trees console readout.
){
  
  registerDoParallel(cores=N_cores)
  if(is.null(New_values)){
    New_values<-0
  }
  
  for(i in 1:length(Drivers)){
    
    if(Use_Gamma){
      print("running with gamma on")
      dir.create(paste(Result_Dest,"/",Which_Parameter,"/",Drivers[i],sep=""),recursive = T)
      foreach(j=Itter)%dopar%{
        for(k in New_values){
          Param_run<-ChangeParameter(Parameter = "useLeafGamma",1,InParamFile = Init_param) %>% 
            ChangeParameter(Parameter=Which_Parameter,NewValue = k, InParamFile = .)
          In_file_run<-ChangeInFile(driver=paste(Driver_loc,"/",Drivers[i],".txt",sep=""),
                                    OutputPrefix = paste(Result_Dest,"/",Which_Parameter,"/",Drivers[i],"/",Which_Parameter,k,"itter",j,sep=""),
                                    ParamFile = Param_run,
                                    InFileLoc = Init_Infile,
                                    param_mod_loc = param_mod,
                                    covfile_loc = covfile
          )
          system(paste(TREES_loc,"<", In_file_run),ignore.stdout=!show_TREES_cout)
          system(paste("echo completed:",Drivers[i],Which_Parameter,k,"itter",j,sep=""))
          unlink(c(Param_run,In_file_run))
        }
      }
    }else{
      print("running with new parameters")
      dir.create(paste(Result_Dest,"/",Which_Parameter,"/",Drivers[i],sep=""),recursive = T)
      print(paste("Created Directory at:",Result_Dest,"/",Which_Parameter,"/",Drivers[i],sep=""))
      foreach(k=New_values)%dopar%{
        Param_run<-ChangeParameter(Parameter=Which_Parameter,NewValue = k, InParamFile = Init_param)
        In_file_run<-ChangeInFile(driver=paste(Driver_loc,"/",Drivers[i],".txt",sep=""),
                                  OutputPrefix = paste(Result_Dest,"/",Which_Parameter,"/",Drivers[i],"/",Which_Parameter,k,"itter",0,sep=""),
                                  ParamFile = Param_run,
                                  InFileLoc = Init_Infile,
                                  param_mod_loc = param_mod,
                                  covfile_loc = covfile
                                  )
        system(paste(TREES_loc,"<", In_file_run),ignore.stdout=!show_TREES_cout)
        system(paste("echo completed:",Drivers[i],Which_Parameter,k,"itter",0,sep=""))
        unlink(c(Param_run,In_file_run))
      }
    }
    
  }
  
  
  
}


###############################################################################################################











