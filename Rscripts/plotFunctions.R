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
                as.expression(bquote("Fine root C (gC "~m^-2 ~")" )),
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
  ylow0 <- min(c(rc0,rc1,rc2,rc3,rc4))
  yhigh0 <- max(c(rc0,rc1,rc2,rc3,rc4))
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  layout(matrix(1:3, ncol = 1, nrow=3), widths = rep(1, 3), heights = c(1.2,1.2,1.2), respect = FALSE)
  par(oma=c(6.1,1.1,2,.1))
  par(mar = c(0.0, 4.8, 0.0, 2.1))
  plot(data$DAY+dayOffset,rc0,type="l",lwd=2, lty=1, col="black", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[4],xaxt="n")
  mtext(title,side=3,cex=1.0, adj=0, line=0.25)
  lines(data$DAY+dayOffset,rc1, lwd=2, lty=1, col="darkgreen")
  lines(data$DAY+dayOffset,rc2, lwd=2, lty=1, col="green4")
  lines(data$DAY+dayOffset,rc3, lwd=2, lty=1, col="green3")
  lines(data$DAY+dayOffset,rc4, lwd=2, lty=1, col="green")
  box(lwd=1.5)
  
  ylow0 <- min(data$YMD..MPa.)
  yhigh0 <- max(data$YPD..MPa.)
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(data$DAY+dayOffset,data$YPD..MPa.,type="l",lwd=2, lty=1, col="black", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[10],xaxt='n')
  lines(data$DAY+dayOffset,data$YPD..MPa., lwd=2, lty=1, col="darkblue")
  lines(data$DAY+dayOffset,data$YMD..MPa., lwd=2, lty=1, col="blue")
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
plotLeafGrowth <- function(subfolder, fname, title, plotname, dayOffset)
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
                as.expression(bquote("SLA ("~m^2 ~kgC^-1~")" )))
  
  leaf<-read.table(paste(subfolder,fname,".leaf",sep=""),header=TRUE,fill=TRUE)
  sim<-read.table(paste(subfolder,fname,".sim",sep=""),header=TRUE)
  
  nrows <- length(sim[,1])
  for(i in 1:nrows)
  {
    sim[i,1] <- dayOffset + (i-1)/48
  }
  #
  #Plot LAI, midday water potential, and GPP for 13 genotypes of maize
  #
  pdf(plotname, width = 5, height= 6.75, useDingbats = F)
  
  al1 <- as.matrix(leaf$Area_Leaf_1)
  al2 <- as.matrix(leaf$Area_Leaf_2)
  al3 <- as.matrix(leaf$Area_Leaf_3)
  al4 <- as.matrix(leaf$Area_Leaf_4)
  al5 <- as.matrix(leaf$Area_Leaf_5)
  al6 <- as.matrix(leaf$Area_Leaf_6)
  al7 <- as.matrix(leaf$Area_Leaf_7)
  ylow0 <- min(c(al1,al2,al3,al4,al5,al6,al7), na.rm=TRUE)
  yhigh0 <- max(c(al1,al2,al3,al4,al5,al6,al7), na.rm=TRUE)
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  layout(matrix(1:3, ncol = 1, nrow=3), widths = rep(1, 3), heights = c(1.2,1.2,1.2), respect = FALSE)
  par(oma=c(6.1,1.1,2,.1))
  par(mar = c(0.0, 4.8, 0.0, 2.1))
  plot(sim[,1],al1,type="l",lwd=2, lty=1, col="black", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[14],xaxt="n")
  mtext(title,side=3,cex=1.0, adj=0, line=0.25)
  lines(sim[,1],al2, lwd=2, lty=1, col="darkgreen")
  lines(sim[,1],al3, lwd=2, lty=1, col="green4")
  lines(sim[,1],al4, lwd=2, lty=1, col="green3")
  lines(sim[,1],al5, lwd=2, lty=1, col="green")
  lines(sim[,1],al6, lwd=2, lty=1, col="green")
  lines(sim[,1],al7, lwd=2, lty=1, col="green")
  box(lwd=1.5)
  
  ylow0 <- min(sim$LAI)
  yhigh0 <- max(sim$LAI)
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(sim[,1],sim$LAI,type="l",lwd=2, lty=1, col="black", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[5],xaxt='n')
  #lines(data$DAY+dayOffset,data$YPD..MPa., lwd=2, lty=1, col="darkblue")
  #lines(data$DAY+dayOffset,data$YMD..MPa., lwd=2, lty=1, col="blue")
  box(lwd=1.5)
  
  ylow0 <- min(sim$SLA)
  yhigh0 <- max(sim$SLA)
  ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
  yhigh <- yhigh0+0.05*(max(abs(ylow0),abs(yhigh0)))
  Axis(side=1, labels=FALSE)
  plot(sim[,1],sim$SLA,type="l",lwd=1, lty=2, col="darkblue", cex.axis=1.5, cex.lab=1.5, 
       ylim=c(ylow,yhigh), ylab=exp_list[15])
  #lines(data$DAY+dayOffset,data$ECRIT..mmol.s.1., lwd=2, lty=1, col="darkred")
  #lines(data$DAY+dayOffset,data$EC..mmol.s.1., lwd=2, lty=1, col="darkblue")
  mtext("Day", side=1, line=3, cex=1.25)
  box(lwd=1.5)
  
  dev.off()
}
