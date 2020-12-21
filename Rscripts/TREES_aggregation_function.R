#
#TREES workshop plotting routines
#  This script is used first to convert half-hourly simulation output to daily
#
#ti	simET	WPlant_K	Soil_Psi	Leaf_Psi	Psi_Crit	Ecrit	Ec	
#RhizFlux0	RhizFlux1	RhizFlux2	RhizFlux3	RhizFlux4	
#Gs	LAI	SLA liveLAI	Rmaint	Rgrowth	leafNSC	stemNSC	rootNSC	
#chloroStarch	chloroSugar waterStress	litterH2O	
#theta0	theta1	theta2	theta3	theta4	thetaRoot	
#Can_Evap	Snowpack	SnowEdef	Vcmax25	Vcmax_sun	Vcmax_shd	
#Jmax25	J_sun	J_shd	Asun	Ashd	Lsun	Lshd	Tsun	Tshd	
#Dsun	Dshd	Ci_sun	Ci_shd	PARsun	PARshd	gs_sun	gs_shd	
#NEE	NPP	R_total	R_ag	R_bg	Rd_sun	Rd_shd	Csapwood	
#FibRootC0	FibRootC1	FibRootC2	FibRootC3	FibRootC4	
#FineRootC0 FineRootC1	FineRootC2	FineRootC3	FineRootC4	
#TotRootC0	TotRootC1	TotRootC2	TotRootC3	TotRootC4	
#FineRootCN0	FineRootCN1	FineRootCNFineRootCN3	FineRootCN4	
#LeafCN	humusC0	humusC1	humusC2	humusC3	humusC4	
#RhizCl0	RhizCl1	RhizCl2	RhizCl3	RhizCl4	
#RhizNl0	RhizNl1	RhizNl2	RhizNl3	RhizNl4	
#AAexudateC0	AAexudateC1	AAexudateC2	AAexudateC3	AAexudateC4	
#SugarExudateC0	SugarExudateC1	SugarExudateC2	SugarExudateC3	SugarExudateC4	
#MicrobC0	MicrobC1	MicrobC2	MicrobC3	MicrobC4	
#MicrobN0	MicrobN1	MicrobN2	MicrobN3	MicrobN4	
#RhizN-	RhizN+	PlantN	PlantNstat	RL ar0	ar1	ar2	ar3	ar4	
#
computeDaily = function(subfolder, fname, simulation, drivers, Ksat)
{
  simulation$WPlant_K <- as.matrix((simulation$WPlant_K))
  simulation[simulation$WPlant_K >Ksat, "WPlant_K"] <- Ksat
  simulation$PLC <- as.matrix((100*(1-simulation$WPlant_K/Ksat)))
  simulation$L <- as.matrix((simulation$Lsun)+(simulation$Lshd))
  simulation$D1 <- as.matrix((simulation$Lsun)*(simulation$Dsun))
  simulation$D2 <- as.matrix((simulation$Lshd)*(simulation$Dshd))
  simulation$DCAN <- as.matrix(((simulation$D1)+(simulation$D2))/(simulation$L))
  simulation$YSOIL <- as.matrix((simulation$Soil_Psi))
  simulation$YPD <- as.matrix((simulation$Leaf_Psi))
  simulation$YMD <- as.matrix((simulation$Leaf_Psi))
  simulation$ECRIT <- as.matrix((simulation$Ecrit)*(simulation$L))
  simulation$NSC <- as.matrix((simulation$leafNSC)+(simulation$stemNSC)+(simulation$rootNSC)
                              +(simulation$chloroStarch)+(simulation$chloroSugar))
  simulation$NSCl <- as.matrix(((simulation$leafNSC)+(simulation$chloroStarch)+(simulation$chloroSugar))/
                                 (max((simulation$liveLAI),(simulation$LAI))/simulation$SLA*10000))
  simulation$NSCs <- as.matrix((simulation$stemNSC)/(simulation$Csapwood))
  simulation$GS <- as.matrix((simulation$Gs))
  simulation$EC <- as.matrix((simulation$Ec)*(simulation$L))
  simulation$GPP <- as.matrix((simulation$Asun)*(simulation$Lsun)+(simulation$Ashd)*(simulation$Lshd))
  simulation$SM01 <- as.matrix(((simulation$theta0)*(simulation$ar0)+(simulation$theta1)*(simulation$ar1))/
                                 ((simulation$ar0)+(simulation$ar1)))
  simulation$SM23 <- as.matrix(((simulation$theta2)*(simulation$ar2)+(simulation$theta3)*(simulation$ar3))/
                                ((simulation$ar2)+(simulation$ar3)))
  simulation$SM4 <- as.matrix(((simulation$theta4)*(simulation$ar4))/((simulation$ar4)))
  simulation$RF01 <- as.matrix(((simulation$RhizFlux0)+(simulation$RhizFlux1))*(simulation$L))
  simulation$RF23 <- as.matrix(((simulation$RhizFlux2)+(simulation$RhizFlux3))*(simulation$L))
  simulation$RF4 <- as.matrix(((simulation$RhizFlux4))*(simulation$L))
  simulation$RF0 <- as.matrix(((simulation$RhizFlux0))*(simulation$L))
  simulation$RF1 <- as.matrix(((simulation$RhizFlux1))*(simulation$L))
  simulation$RF2 <- as.matrix(((simulation$RhizFlux2))*(simulation$L))
  simulation$RF3 <- as.matrix(((simulation$RhizFlux3))*(simulation$L))
  simulation$RF4 <- as.matrix(((simulation$RhizFlux4))*(simulation$L))
  simulation$FineRC01 <- as.matrix((simulation$FibRootC0)+(simulation$FibRootC1)+(simulation$FineRootC0)+(simulation$FineRootC1))
  simulation$FineRC23 <- as.matrix((simulation$FibRootC2)+(simulation$FineRootC2)+(simulation$FibRootC3)+(simulation$FineRootC3))
  simulation$FineRC4 <- as.matrix((simulation$FibRootC4)+(simulation$FineRootC4))
  simulation$FineRC0 <- as.matrix((simulation$FibRootC0)+(simulation$FineRootC0))
  simulation$FineRC1 <- as.matrix((simulation$FibRootC1)+(simulation$FineRootC1))
  simulation$FineRC2 <- as.matrix((simulation$FibRootC2)+(simulation$FineRootC2))
  simulation$FineRC3 <- as.matrix((simulation$FibRootC3)+(simulation$FineRootC3))
  simulation$FineRC4 <- as.matrix((simulation$FibRootC4)+(simulation$FineRootC4))
  simulation$TotalRC01 <- as.matrix((simulation$TotRootC0)+(simulation$TotRootC1))
  simulation$TotalRC23 <- as.matrix((simulation$TotRootC2)+(simulation$TotRootC3))
  simulation$TotalRC4 <- as.matrix((simulation$TotRootC4))
  simulation$FineRootC <- as.matrix((simulation$FineRC01)+(simulation$FineRC23)+(simulation$FineRC4))
  simulation$TotalRootC <- as.matrix((simulation$TotalRC01)+(simulation$TotalRC23)+(simulation$TotalRC4))
  simulation$RootC <- as.matrix((simulation$TotalRootC))
  simulation$NSCr <- as.matrix((simulation$rootNSC)/(simulation$RootC))
  simulation$MicrobN01 <- as.matrix((simulation$MicrobN0)+(simulation$MicrobN1))
  simulation$MicrobN23 <- as.matrix((simulation$MicrobN2)+(simulation$MicrobN3))
  simulation$MicrobN4 <- as.matrix((simulation$MicrobN4))
  simulation$RhizNl01 <- as.matrix((simulation$RhizNl0)+(simulation$RhizNl1))
  simulation$RhizNl23 <- as.matrix((simulation$RhizNl2)+(simulation$RhizNl3))
  simulation$RhizNl4 <- as.matrix((simulation$RhizNl4))
  simulation$AR01 <- as.matrix((simulation$ar0)+(simulation$ar1))
  simulation$AR23 <- as.matrix((simulation$ar2)+(simulation$ar3))
  simulation$AR4 <- as.matrix((simulation$ar4))
  simulation$reprod <- as.matrix((simulation$reproduction))
  
  DCAN_conversion <- 1.0/4.0
  PLC_conversion <- 1.0/4.0
  YSOIL_conversion <- 1.0/2.0
  YPD_conversion <- 1.0/2.0
  YMD_conversion <- 1.0/2.0
  ECRIT_conversion <- 1.0/2.0
  leafNSC_conversion <- 1.0/4.0/10.0
  stemNSC_conversion <- 1.0/4.0/10.0
  rootNSC_conversion <- 1.0/4.0/10.0
  NSCr_conversion <- 1.0/4.0*100.0
  GS_conversion <- 1.0/2.0*1000
  EC_conversion <- 1.0/2.0
  GPP_conversion <- 1.0/48.0
  SM01_conversion <- 1.0/48.0
  SM23_conversion <- 1.0/48.0
  SM4_conversion <- 1.0/48.0
  theta0_conversion <- 1.0/48.0
  theta1_conversion <- 1.0/48.0
  theta2_conversion <- 1.0/48.0
  theta3_conversion <- 1.0/48.0
  theta4_conversion <- 1.0/48.0
  theta5_conversion <- 1.0/48.0
  RF01_conversion <- 1.0/2.0
  RF23_conversion <- 1.0/2.0
  RF4_conversion <- 1.0/2.0
  RF0_conversion <- 1.0/2.0
  RF1_conversion <- 1.0/2.0
  RF2_conversion <- 1.0/2.0
  RF3_conversion <- 1.0/2.0
  RF4_conversion <- 1.0/2.0
  FineRC01_conversion <- 1.0/48.0/10.0
  FineRC23_conversion <- 1.0/48.0/10.0
  FineRC4_conversion <- 1.0/48.0/10.0
  FineRC0_conversion <- 1.0/48.0/10.0
  FineRC1_conversion <- 1.0/48.0/10.0
  FineRC2_conversion <- 1.0/48.0/10.0
  FineRC3_conversion <- 1.0/48.0/10.0
  FineRC4_conversion <- 1.0/48.0/10.0
  MicrobN01_conversion <- 1.0/48.0/10.0
  MicrobN23_conversion <- 1.0/48.0/10.0
  MicrobN4_conversion <- 1.0/48.0/10.0
  MicrobN0_conversion <- 1.0/48.0/10.0
  MicrobN1_conversion <- 1.0/48.0/10.0
  MicrobN2_conversion <- 1.0/48.0/10.0
  MicrobN3_conversion <- 1.0/48.0/10.0
  MicrobN4_conversion <- 1.0/48.0/10.0
  RhizNl01_conversion <- 1.0/48.0/10.0
  RhizNl23_conversion <- 1.0/48.0/10.0
  RhizNl4_conversion <- 1.0/48.0/10.0
  PlantN_conversion <- 1.0/48.0/10.0
  LAI_conversion <- 1.0/48.0
  SLA_conversion <- 1.0/48.0
  RLA_conversion <- 1.0/48.0
  AR01_conversion <- 1.0/48
  AR23_conversion <- 1.0/48.0
  AR4_conversion <- 1.0/48.0
  ar0_conversion <- 1.0/48
  ar1_conversion <- 1.0/48
  ar2_conversion <- 1.0/48
  ar3_conversion <- 1.0/48
  ar4_conversion <- 1.0/48
  reprod_conversion <- 1.0/48.0/10.0
  
  flux_list<-c("DCAN","PLC","YSOIL","YPD","YMD","ECRIT","leafNSC","stemNSC","rootNSC","GS","EC","GPP",
               "theta0","theta1","theta2","theta3","theta4","RF0","RF1","RF2","RF3","RF4",
               "FineRC0","FineRC1","FineRC2","FineRC3","FineRC4",
               "MicrobN0","MicrobN1","MicrobN2","MicrobN3","MicrobN4",
               "PlantN","LAI","SLA","RLA","ar0","ar1","ar2","ar3","ar4","reprod")
  
#start and end times
  daylow <- c(12,12,3.5,3.5,13,13,12,12,12,13,13,0,
              0,0,0,0,0,13,13,13,13,13,
              0,0,0,0,0,
              0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0)
  dayhigh <- c(13.5,13.5,4,4,13.5,13.5,13.5,13.5,13.5,13.5,13.5,23.5,
               23.5,23.5,23.5,23.5,23.5,13.5,13.5,13.5,13.5,13.5,
               23.5,23.5,23.5,23.5,23.5,
               23.5,23.5,23.5,23.5,23.5,
               23.5,23.5,23.5,23.5,23.5,23.5,23.5,23.5,23.5,23.5)
  
#allocate array for output
  nrows <- length(simulation[,1])/48
  result_array <- array(data=0,dim=c(nrows,43))
  for(i in 1:nrows)
  {
    result_array[i,1] <- i
  }
  
  for(j in 1:length(flux_list))
  {
#Binding together the driver file and simulation columns
    flux_bound<-cbind(drivers[,1:2],simulation[,flux_list[j]])
    colnames(flux_bound)<-c("jday","Hour","simulated")
    midday_flux<-subset(flux_bound,flux_bound$Hour>=daylow[j] & flux_bound$Hour <=dayhigh[j])
#Aggregating to a daily timestep
    flux_agg<-aggregate(midday_flux,by=list(midday_flux$jday),FUN=sum,na.rm=TRUE)
#Making the proper conversion 
    flux_agg$simulated<-flux_agg$simulated*get(paste(flux_list[j],"_conversion",sep=""))
    result_array[,j+1]<-as.vector(flux_agg$simulated)
  }
  
  colnames(result_array)<-c("DAY","DCAN (kPa)","PLC","YSOIL (MPa)","YPD (MPa)", "YMD (MPa)", 
                            "ECRIT (mmol s-1)","Leaf NSC (gC m-2)","Stem NSC (gC m-2)","Root NSC (gC m-2)",
                            "GS (mmol m-2 s-1)","EC (mmol s-1)","GPP (mmol m-2 s-1)",
                            "theta0","theta1","theta2","theta3","theta4",
                            "RHIZF0 (mmol s-1)","RHIZF1 (mmol s-1)","RHIZF2 (mmol s-1)","RHIZF3 (mmol s-1)","RHIZF4 (mmol s-1)",
                            "FinRootC0 (gC m-2 grd)","FinRootC1 (gC m-2 grd)","FinRootC2 (gC m-2 grd)","FinRootC3 (gC m-2 grd)","FinRootC4 (gC m-2 grd)",
                            "MicrobN0 (gN m-2 grd)","MicrobN1 (gN m-2 grd)","MicrobN2 (gN m-2 grd)","MicrobN3 (gN m-2 grd)","MicrobN4 (gN m-2 grd)",
                            "PlantN (gN m-2 grd)", "LAI (m2 m-2)", "SLA (m-2 kgC)","RLA (m2 m-2)",
                            "RA0 (m2 m-2","RA1 (m2 m-2)","RA3 (m2 m-2)","RA4 (m2 m-2)","RA4 (m2 m-2)","Reproduction (gC m-2)")
  write.csv(result_array, file=paste(subfolder,fname,"_midday.csv"))
}

