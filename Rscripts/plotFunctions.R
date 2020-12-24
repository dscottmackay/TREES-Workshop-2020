#
#
#
plotDaily <- function(data, title, plotname, dayOffset)
{
  #
  #Define some y-axis labels
  #
  exp_list <- c(as.expression(bquote(Psi[MD]~"(MPa)" )),
                as.expression(bquote(italic(E)[C]~"or"~italic(E)[Crit]~"("~"mmol" ~s^-1~")" )),
                as.expression(bquote("Relative root area" )),
                as.expression(bquote("Leaf and FR C (gC "~m^-2 ~")" )),
                as.expression(bquote("LAI ("~m^2 ~m^-2~")" )),
                as.expression(bquote("Reproduction (gC "~m^-2 ~")" )),
                as.expression(bquote("Plant N (gN "~m^-2 ~")" )),
                as.expression(bquote("NSC (gC "~m^-2 ~")" )),
                as.expression(bquote(italic(E)[C]~"or"~italic(F)[Rhiz]~"("~"mmol" ~s^-1~")" )),
                as.expression(bquote(Psi[Leaf]~"(MPa)" )),
                as.expression(bquote("PLC (%)")),
                as.expression(bquote("GPP (umol "~m^-2 ~"gnd" ~s^-1 ~")" )),
                as.expression(bquote(italic(E)[C] ~" (mmol "~m^-2 ~"gnd" ~s^-1 ~")" )))
  
  #
  #Plot LAI, midday water potential, and GPP for 13 genotypes of maize
  #
  pdf(plotname, width = 5, height= 6.75, useDingbats = F)
  
  rc0 <- as.matrix(data$FinRootC0..gC.m.2.grd.)
  rc1 <- as.matrix(data$FinRootC1..gC.m.2.grd.)
  rc2 <- as.matrix(data$FinRootC2..gC.m.2.grd.)
  rc3 <- as.matrix(data$FinRootC3..gC.m.2.grd.)
  rc4 <- as.matrix(data$FinRootC4..gC.m.2.grd.)
  rc <- as.matrix(rc0+rc1+rc2+rc3+rc4)
  lc <- as.matrix((data$LAI..m2.m.2.)/(data$SLA..m.2.kgC.)*1000)
  ylow0 <- min(c(rc0,rc1,rc2,rc3,rc4,lc))
  yhigh0 <- max(c(rc0,rc1,rc2,rc3,rc4,lc))
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  layout(matrix(1:3, ncol = 1, nrow=3), widths = rep(1, 3), heights = c(1.2,1.2,1.2), respect = FALSE)
  par(oma=c(6.1,1.1,2,.1))
  par(mar = c(0.0, 4.8, 0.0, 2.1))
  plot(data$DAY+dayOffset,rc0,type="l",lwd=2, lty=1, col="black", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[4],xaxt="n")
  mtext(title,side=3,cex=1.0, adj=0, line=0.25)
  lines(data$DAY+dayOffset,rc1, lwd=2, lty=1, col="tan4")
  lines(data$DAY+dayOffset,rc2, lwd=2, lty=1, col="tan3")
  lines(data$DAY+dayOffset,rc3, lwd=2, lty=1, col="tan2")
  lines(data$DAY+dayOffset,rc4, lwd=2, lty=1, col="tan1")
  lines(data$DAY+dayOffset,rc, lwd=2, lty=2, col="tan4")
  lines(data$DAY+dayOffset,lc, lwd=2, lty=1, col="darkgreen")
  box(lwd=1.5)
  
  #ylow0 <- min(data$YMD..MPa.)
  #yhigh0 <- max(data$YPD..MPa.)
  #ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  #yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  #Axis(side=1, labels=FALSE)
  #plot(data$DAY+dayOffset,data$YPD..MPa.,type="l",lwd=2, lty=1, col="black", cex.axis=1.5, cex.lab=1.5, 
  #     ylim=c(ylow,yhigh), ylab=exp_list[10],xaxt='n')
  #lines(data$DAY+dayOffset,data$YPD..MPa., lwd=2, lty=1, col="darkblue")
  #lines(data$DAY+dayOffset,data$YMD..MPa., lwd=2, lty=1, col="blue")
  #box(lwd=1.5)
  
  ylow0 <- min(data$GPP..mmol.m.2.s.1.)
  yhigh0 <- max(data$GPP..mmol.m.2.s.1.)
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(data$DAY+dayOffset,data$GPP..mmol.m.2.s.1.,type="l",lwd=2, lty=1, col="black", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[12],xaxt='n')
  box(lwd=1.5)
  
  ylow0 <- min(data$EC..mmol.s.1.)
  yhigh0 <- max(data$ECRIT..mmol.s.1.)
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(data$DAY+dayOffset,data$EC..mmol.s.1.,type="l",lwd=1, lty=2, col="darkblue", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[13])
  lines(data$DAY+dayOffset,data$ECRIT..mmol.s.1., lwd=2, lty=1, col="darkred")
  lines(data$DAY+dayOffset,data$EC..mmol.s.1., lwd=2, lty=1, col="darkblue")
  mtext("Day", side=1, line=3, cex=1.25)
  box(lwd=1.5)
  
  dev.off()
}

#
#
#
plotDaily2 <- function(data, data2, title, plotname, dayOffset)
{
  #
  #Define some y-axis labels
  #
  exp_list <- c(as.expression(bquote(Psi[MD]~"(MPa)" )),
                as.expression(bquote(italic(E)[C]~"or"~italic(E)[Crit]~"("~"mmol" ~s^-1~")" )),
                as.expression(bquote("Relative root area" )),
                as.expression(bquote("Leaf and FR C (gC "~m^-2 ~")" )),
                as.expression(bquote("LAI ("~m^2 ~m^-2~")" )),
                as.expression(bquote("Reproduction (gC "~m^-2 ~")" )),
                as.expression(bquote("Plant N (gN "~m^-2 ~")" )),
                as.expression(bquote("NSC (gC "~m^-2 ~")" )),
                as.expression(bquote(italic(E)[C]~"or"~italic(F)[Rhiz]~"("~"mmol" ~s^-1~")" )),
                as.expression(bquote(Psi[Leaf]~"(MPa)" )),
                as.expression(bquote("PLC (%)")),
                as.expression(bquote("GPP (umol "~m^-2 ~"gnd" ~s^-1 ~")" )),
                as.expression(bquote(italic(E)[C] ~" (mmol "~m^-2 ~"gnd" ~s^-1 ~")" )))
  
  #
  #Plot LAI, midday water potential, and GPP for 13 genotypes of maize
  #
  pdf(plotname, width = 5, height= 6.75, useDingbats = F)
  
  nscl <- as.matrix(data$Leaf.NSC..gC.m.2.)
  nscs <- as.matrix(data$Stem.NSC..gC.m.2.)
  nscr <- as.matrix(data$Root.NSC..gC.m.2.)
  nscl2 <- as.matrix(data2$Leaf.NSC..gC.m.2.)
  nscs2 <- as.matrix(data2$Stem.NSC..gC.m.2.)
  nscr2 <- as.matrix(data2$Root.NSC..gC.m.2.)
  ylow0 <- min(c(nscl,nscs,nscr,nscl2,nscs2,nscr2))
  yhigh0 <- max(c(nscl,nscs,nscr,nscl2,nscs2,nscr2))
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  layout(matrix(1:3, ncol = 1, nrow=3), widths = rep(1, 3), heights = c(1.2,1.2,1.2), respect = FALSE)
  par(oma=c(6.1,1.1,2,.1))
  par(mar = c(0.0, 4.8, 0.0, 2.1))
  plot(data$DAY+dayOffset,nscl,type="l",lwd=2, lty=1, col="darkgreen", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[8],xaxt="n")
  mtext(title,side=3,cex=1.0, adj=0, line=0.25)
  lines(data$DAY+dayOffset,nscs, lwd=2, lty=1, col="tan4")
  lines(data$DAY+dayOffset,nscr, lwd=2, lty=1, col="tan2")
  lines(data$DAY+dayOffset,nscl2, lwd=2, lty=2, col="darkgreen")
  lines(data$DAY+dayOffset,nscs2, lwd=2, lty=2, col="tan4")
  lines(data$DAY+dayOffset,nscr2, lwd=2, lty=2, col="tan2")
  box(lwd=1.5)
  
  ylow0 <- min(c(data$YMD..MPa.,data2$YMD..MPa))
  yhigh0 <- max(c(data$YPD..MPa.,data2$YMD..MPa))
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(data$DAY+dayOffset,data$YPD..MPa.,type="l",lwd=2, lty=1, col="black", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[10],xaxt='n')
  lines(data$DAY+dayOffset,data$YPD..MPa., lwd=2, lty=1, col="darkblue")
  lines(data$DAY+dayOffset,data$YMD..MPa., lwd=2, lty=1, col="blue")
  lines(data$DAY+dayOffset,data2$YPD..MPa., lwd=2, lty=2, col="darkblue")
  lines(data$DAY+dayOffset,data2$YMD..MPa., lwd=2, lty=2, col="blue")
  box(lwd=1.5)
  
  #ylow0 <- min(c(data$EC..mmol.s.1.,data2$EC..mmol.s.1.))
  #yhigh0 <- max(c(data$ECRIT..mmol.s.1.,data2$ECRIT..mmol.s.1.))
  #ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  #yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  #Axis(side=1, labels=FALSE)
  #plot(data$DAY+dayOffset,data$EC..mmol.s.1.,type="l",lwd=1, lty=2, col="darkblue", cex.axis=1.5, cex.lab=1.5, 
  #     ylim=c(ylow,yhigh), ylab=exp_list[13])
  #lines(data$DAY+dayOffset,data$ECRIT..mmol.s.1., lwd=2, lty=1, col="darkred")
  #lines(data$DAY+dayOffset,data$EC..mmol.s.1., lwd=2, lty=1, col="darkblue")
  #lines(data$DAY+dayOffset,data2$ECRIT..mmol.s.1., lwd=2, lty=2, col="darkred")
  #lines(data$DAY+dayOffset,data2$EC..mmol.s.1., lwd=2, lty=2, col="darkblue")
  #mtext("Day", side=1, line=3, cex=1.25)
  #box(lwd=1.5)
  
  ylow0 <- min(c(data$PLC,data2$PLC))
  yhigh0 <- max(c(data$PLC,data2$PLC))
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(data$DAY+dayOffset,data$PLC,type="l",lwd=2, lty=1, col="darkred", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[11])
  lines(data$DAY+dayOffset,data2$PLC, lwd=2, lty=2, col="darkred")
  mtext("Day", side=1, line=3, cex=1.25)
  box(lwd=1.5)
  
  dev.off()
}

plotDaily3 <- function(data, data2, data3, title, plotname, dayOffset)
{
  #
  #Define some y-axis labels
  #
  exp_list <- c(as.expression(bquote(Psi[MD]~"(MPa)" )),
                as.expression(bquote(italic(E)[C]~"or"~italic(E)[Crit]~"("~"mmol" ~s^-1~")" )),
                as.expression(bquote("Relative root area" )),
                as.expression(bquote("Leaf and FR C (gC "~m^-2 ~")" )),
                as.expression(bquote("LAI ("~m^2 ~m^-2~")" )),
                as.expression(bquote("Reproduction (gC "~m^-2 ~")" )),
                as.expression(bquote("Plant N (gN "~m^-2 ~")" )),
                as.expression(bquote("NSC (gC "~m^-2 ~")" )),
                as.expression(bquote(italic(E)[C]~"or"~italic(F)[Rhiz]~"("~"mmol" ~s^-1~")" )),
                as.expression(bquote(Psi[Leaf]~"(MPa)" )),
                as.expression(bquote("PLC (%)")),
                as.expression(bquote("GPP (umol "~m^-2 ~"gnd" ~s^-1 ~")" )),
                as.expression(bquote(italic(E)[C] ~" (mmol "~m^-2 ~"gnd" ~s^-1 ~")" )))
  
  #
  #Plot LAI, midday water potential, and GPP for 13 genotypes of maize
  #
  pdf(plotname, width = 5, height= 6.75, useDingbats = F)
  
  rc0 <- as.matrix(data$FinRootC0..gC.m.2.grd.)
  rc1 <- as.matrix(data$FinRootC1..gC.m.2.grd.)
  rc2 <- as.matrix(data$FinRootC2..gC.m.2.grd.)
  rc3 <- as.matrix(data$FinRootC3..gC.m.2.grd.)
  rc4 <- as.matrix(data$FinRootC4..gC.m.2.grd.)
  rootc <- as.matrix(rc0+rc1+rc2+rc3+rc4)
  leafc <- as.matrix((data$LAI..m2.m.2.)/(data$SLA..m.2.kgC.)*1000)
  rc0 <- as.matrix(data2$FinRootC0..gC.m.2.grd.)
  rc1 <- as.matrix(data2$FinRootC1..gC.m.2.grd.)
  rc2 <- as.matrix(data2$FinRootC2..gC.m.2.grd.)
  rc3 <- as.matrix(data2$FinRootC3..gC.m.2.grd.)
  rc4 <- as.matrix(data2$FinRootC4..gC.m.2.grd.)
  rootc2 <- as.matrix(rc0+rc1+rc2+rc3+rc4)
  leafc2<- as.matrix((data2$LAI..m2.m.2.)/(data2$SLA..m.2.kgC.)*1000)
  rc0 <- as.matrix(data3$FinRootC0..gC.m.2.grd.)
  rc1 <- as.matrix(data3$FinRootC1..gC.m.2.grd.)
  rc2 <- as.matrix(data3$FinRootC2..gC.m.2.grd.)
  rc3 <- as.matrix(data3$FinRootC3..gC.m.2.grd.)
  rc4 <- as.matrix(data3$FinRootC4..gC.m.2.grd.)
  rootc3 <- as.matrix(rc0+rc1+rc2+rc3+rc4)
  leafc3<- as.matrix((data3$LAI..m2.m.2.)/(data3$SLA..m.2.kgC.)*1000)
  ylow0 <- min(c(rootc,leafc,rootc2,leafc2,rootc3,leafc3))
  yhigh0 <- max(c(rootc,leafc,rootc2,leafc2,rootc3,leafc3))
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  layout(matrix(1:3, ncol = 1, nrow=3), widths = rep(1, 3), heights = c(1.2,1.2,1.2), respect = FALSE)
  par(oma=c(6.1,1.1,2,.1))
  par(mar = c(0.0, 4.8, 0.0, 2.1))
  plot(data$DAY+dayOffset,leafc,type="l",lwd=2, lty=1, col="darkgreen", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[4],xaxt="n")
  mtext(title,side=3,cex=1.0, adj=0, line=0.25)
  lines(data$DAY+dayOffset,leafc2, lwd=2, lty=2, col="darkgreen")
  lines(data$DAY+dayOffset,leafc3, lwd=2, lty=5, col="darkgreen")
  lines(data$DAY+dayOffset,rootc, lwd=2, lty=1, col="tan4")
  lines(data$DAY+dayOffset,rootc2, lwd=2, lty=2, col="tan4")
  lines(data$DAY+dayOffset,rootc3, lwd=2, lty=5, col="tan4")
  #plot C11 
  points(c(186,205,227,247),c(2.159,4.513,4.201,4.706)/33*1000,
         pch=1,col="darkgreen",cex=1.25,lwd=1.25)
  #plot C22
  points(c(186,205,227,247),c(2.007,3.83,NA,4.607)/33*1000,
         pch=1,col="darkgreen",cex=1.25,lwd=1.25)
  #plot D31
  points(c(186,205,227,247),c(NA,3.307,4.446,4.443)/33*1000,
         pch=1,col="darkgreen",cex=1.25,lwd=1.25)
  #plot D43
  points(c(186,205,227,247),c(NA,3.551,NA,4.357)/33*1000,
         pch=1,col="darkgreen",cex=1.25,lwd=1.25)
  legend("bottomright",lty=c(1,2,5),lwd=c(2,2,2), cex=0.7,
         c("B73 - Temperate origin","MO18W - Mixed origin","CML103 - Tropical origin"))
  box(lwd=1.5)
  
  ylow0 <- min(c(data$GPP..mmol.m.2.s.1.,data2$GPP..mmol.m.2.s.1.,data3$GPP..mmol.m.2.s.1.))
  yhigh0 <- max(c(data$GPP..mmol.m.2.s.1.,data2$GPP..mmol.m.2.s.1.,data3$GPP..mmol.m.2.s.1.))
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(data$DAY+dayOffset,data$GPP..mmol.m.2.s.1.,type="l",lwd=2, lty=1, col="black", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[12],xaxt='n')
  lines(data$DAY+dayOffset,data2$GPP..mmol.m.2.s.1., lwd=2, lty=2, col="black")
  lines(data$DAY+dayOffset,data3$GPP..mmol.m.2.s.1., lwd=2, lty=5, col="black")
  box(lwd=1.5)
  
  #ylow0 <- min(c(data$EC..mmol.s.1.,data2$EC..mmol.s.1.,data3$EC..mmol.s.1.))
  #yhigh0 <- max(c(data$ECRIT..mmol.s.1.,data2$ECRIT..mmol.s.1.,data3$EC..mmol.s.1.))
  #ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  #yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  #Axis(side=1, labels=FALSE)
  #plot(data$DAY+dayOffset,data$EC..mmol.s.1.,type="l",lwd=2, lty=1, col="darkblue", cex.axis=1.5, cex.lab=1.5, 
  #     ylim=c(ylow,yhigh), ylab=exp_list[13])
  #lines(data$DAY+dayOffset,data$ECRIT..mmol.s.1., lwd=2, lty=1, col="darkred")
  #lines(data$DAY+dayOffset,data2$EC..mmol.s.1., lwd=2, lty=2, col="darkblue")
  #lines(data$DAY+dayOffset,data2$ECRIT..mmol.s.1., lwd=2, lty=2, col="darkred")
  #lines(data$DAY+dayOffset,data3$EC..mmol.s.1., lwd=2, lty=5, col="darkblue")
  #lines(data$DAY+dayOffset,data3$ECRIT..mmol.s.1., lwd=2, lty=5, col="darkred")
  #mtext("Day", side=1, line=3, cex=1.25)
  #box(lwd=1.5)
  
  ylow0 <- min(c(data$PLC,data2$PLC,data3$PLC))
  yhigh0 <- max(c(data$PLC,data2$PLC,data3$PLC))
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(data$DAY+dayOffset,data$PLC,type="l",lwd=2, lty=1, col="darkred", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[11])
  lines(data$DAY+dayOffset,data2$PLC, lwd=2, lty=2, col="darkred")
  lines(data$DAY+dayOffset,data3$PLC, lwd=2, lty=5, col="darkred")
  mtext("Day", side=1, line=3, cex=1.25)
  box(lwd=1.5)
  
  dev.off()
}

#
#
#
plotLeafGrowth <- function(subfolder, fname1, fname2, title, plotname, dayOffset)
{
  #
  #Define some y-axis labels
  #
  exp_list <- c(as.expression(bquote(Psi[MD]~"(MPa)" )),
                as.expression(bquote(italic(E)[C]~"or"~italic(E)[Crit]~"("~"mmol" ~s^-1~")" )),
                as.expression(bquote("Relative root area" )),
                as.expression(bquote("Fine root C (gC "~m^-2 ~")" )),
                as.expression(bquote("LAI ("~m^2 ~m^-2~")" )),
                as.expression(bquote("Reproduction (gC "~m^-2 ~")" )),
                as.expression(bquote("Plant N (gN "~m^-2 ~")" )),
                as.expression(bquote("NSC (gC "~m^-2 ~")" )),
                as.expression(bquote(italic(E)[C]~"or"~italic(F)[Rhiz]~"("~"mmol" ~s^-1~")" )),
                as.expression(bquote(Psi[Leaf]~"(MPa)" )),
                as.expression(bquote("PLC (%)")),
                as.expression(bquote("GPP (umol "~m^-2 ~"gnd" ~s^-1 ~")" )),
                as.expression(bquote(italic(E)[C] ~" (mmol "~m^-2 ~"gnd" ~s^-1 ~")" )),
                as.expression(bquote("Leaf Area ("~cm^2~")" )),
                as.expression(bquote("SLA ("~m^2 ~kgC^-1~")" )),
                as.expression(bquote("C:N ("~kgC^2 ~kgN^-1~")" )))
  
  leaf1<-read.table(paste(subfolder,fname1,".leaf",sep=""),header=TRUE,fill=TRUE)
  leaf2<-read.table(paste(subfolder,fname2,".leaf",sep=""),header=TRUE,fill=TRUE)
  sim1<-read.table(paste(subfolder,fname1,".sim",sep=""),header=TRUE)
  sim2<-read.table(paste(subfolder,fname2,".sim",sep=""),header=TRUE)
  
  nrows <- length(sim1[,1])
  for(i in 1:nrows)
  {
    sim1[i,1] <- dayOffset + (i-1)/48
  }
  #
  #Plot LAI, midday water potential, and GPP for 13 genotypes of maize
  #
  pdf(plotname, width = 5, height= 6.75, useDingbats = F)
  
  al1 <- as.matrix(leaf1$Area_Leaf_1)
  al2 <- as.matrix(leaf1$Area_Leaf_2)
  al3 <- as.matrix(leaf1$Area_Leaf_3)
  al4 <- as.matrix(leaf1$Area_Leaf_4)
  al5 <- as.matrix(leaf1$Area_Leaf_5)
  al6 <- as.matrix(leaf1$Area_Leaf_6)
  al7 <- as.matrix(leaf1$Area_Leaf_7)
  all1 <- as.matrix(leaf2$Area_Leaf_1)
  all2 <- as.matrix(leaf2$Area_Leaf_2)
  all3 <- as.matrix(leaf2$Area_Leaf_3)
  all4 <- as.matrix(leaf2$Area_Leaf_4)
  all5 <- as.matrix(leaf2$Area_Leaf_5)
  all6 <- as.matrix(leaf2$Area_Leaf_6)
  all7 <- as.matrix(leaf2$Area_Leaf_7)
  ylow0 <- min(c(al1,al2,al3,al4,al5,al6,al7,all1,all2,all3,all4,all5,all6,all7), na.rm=TRUE)
  yhigh0 <- max(c(al1,al2,al3,al4,al5,al6,al7,all1,all2,all3,all4,all5,all6,all7), na.rm=TRUE)
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  layout(matrix(1:3, ncol = 1, nrow=3), widths = rep(1, 3), heights = c(1.2,1.2,1.2), respect = FALSE)
  par(oma=c(6.1,1.1,2,.1))
  par(mar = c(0.0, 4.8, 0.0, 2.1))
  plot(sim1[,1],al1,type="l",lwd=2, lty=1, col="darkgreen", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[14],xaxt="n")
  mtext(title,side=3,cex=1.0, adj=0, line=0.25)
  lines(sim1[,1],al2, lwd=2, lty=1, col="darkgreen")
  lines(sim1[,1],al3, lwd=2, lty=1, col="darkgreen")
  lines(sim1[,1],al4, lwd=2, lty=1, col="darkgreen")
  lines(sim1[,1],al5, lwd=2, lty=1, col="darkgreen")
  lines(sim1[,1],al6, lwd=2, lty=1, col="darkgreen")
  lines(sim1[,1],al7, lwd=2, lty=1, col="darkgreen")
  lines(sim1[,1],all1, lwd=2, lty=1, col="green")
  lines(sim1[,1],all2, lwd=2, lty=1, col="green")
  lines(sim1[,1],all3, lwd=2, lty=1, col="green")
  lines(sim1[,1],all4, lwd=2, lty=1, col="green")
  lines(sim1[,1],all5, lwd=2, lty=1, col="green")
  lines(sim1[,1],all6, lwd=2, lty=1, col="green")
  lines(sim1[,1],all7, lwd=2, lty=1, col="green")
  box(lwd=1.5)
  
  al1 <- as.matrix(sim1$FineRootCN0)
  al2 <- as.matrix(sim1$FineRootCN1)
  al3 <- as.matrix(sim1$FineRootCN2)
  al4 <- as.matrix(sim1$FineRootCN3)
  al5 <- as.matrix(sim1$FineRootCN4)
  al6 <- as.matrix(sim1$LeafCN)
  all1 <- as.matrix(sim2$FineRootCN0)
  all2 <- as.matrix(sim2$FineRootCN1)
  all3 <- as.matrix(sim2$FineRootCN2)
  all4 <- as.matrix(sim2$FineRootCN3)
  all5 <- as.matrix(sim2$FineRootCN4)
  all6 <- as.matrix(sim2$LeafCN)
  ylow0 <- min(c(al1,al2,al3,al4,al5,al6,all1,all2,all3,all4,all5,all6), na.rm=TRUE)
  yhigh0 <- max(c(al1,al2,al3,al4,al5,al6,all1,all2,all3,all4,all5,all6), na.rm=TRUE)
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(sim1[,1],al1,type="l",lwd=2, lty=1, col="tan4", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[16],xaxt='n')
  lines(sim1[,1],al2, lwd=2, lty=1, col="tan4")
  lines(sim1[,1],al3, lwd=2, lty=1, col="tan4")
  lines(sim1[,1],al4, lwd=2, lty=1, col="tan4")
  lines(sim1[,1],al5, lwd=2, lty=1, col="tan4")
  lines(sim1[,1],al6, lwd=2, lty=1, col="darkgreen")
  lines(sim1[,1],all1, lwd=2, lty=1, col="tan2")
  lines(sim1[,1],all2, lwd=2, lty=1, col="tan2")
  lines(sim1[,1],all3, lwd=2, lty=1, col="tan2")
  lines(sim1[,1],all4, lwd=2, lty=1, col="tan2")
  lines(sim1[,1],all5, lwd=2, lty=1, col="tan2")
  lines(sim1[,1],all6, lwd=2, lty=1, col="green")
  #lines(data$DAY+dayOffset,data$YMD..MPa., lwd=2, lty=1, col="blue")
  box(lwd=1.5)
  
  ylow0 <- min(c(sim1$SLA,sim2$SLA))
  yhigh0 <- max(c(sim1$SLA,sim2$SLA))
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(sim1[,1],sim1$SLA,type="l",lwd=2, lty=1, col="darkgreen", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[15])
  lines(sim1[,1],sim2$SLA, lwd=2, lty=1, col="green")
  #lines(data$DAY+dayOffset,data$EC..mmol.s.1., lwd=2, lty=1, col="darkblue")
  mtext("Day", side=1, line=3, cex=1.25)
  box(lwd=1.5)
  
  dev.off()
}
