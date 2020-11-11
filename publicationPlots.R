#
# This file uses the output of hiv.R to produce figure 4 (and SI equivalents)
#
filename<-"blackmale"
#filename<-"blackfemale"
#filename<-"whitemale"
#filename<-"whitefemale"
load(paste0(filename,".Rdata"))
source('hivFunctions.R')
yearbyyear<-yearbyyear.full
uprq<-0.95 # set upper and lower quantiles for 90% CI
lwrq<-0.05
dis<-c("HCV","HIV","HSV")
birth.years<-c(paste0("-",byrs[2]-1),byrs[2:ags])
yrs<-2002.5+1:(tps-1)
sps<-(iters-burnin)/thinning
#
# Open plot file
#
pdf(paste0("fig_",filename,".pdf"),paper="a4",width=9.5,height=12,pointsize=10)
par(mfrow=c(4,3))
#
# Model fit of the 3 diseases
#
a.col<-c("#00990055","#FF000055","#0000FF55")
cols<-c("#009900","#FF0000","#0000FF")
disease.status<-c("HCV","HIV","HSV2")
joins<-list(5:8,c(3,4,7,8),2*1:4)
data.years<-c(1,5,11)
for (b in 1:3) {
  mid<-matrix(NA,ags,length(joins))
  lwr<-matrix(NA,ags,length(joins))
  upr<-matrix(NA,ags,length(joins))
  for (a in 1:length(joins)) {
    mid[,a]<-apply(apply(x.stored[data.years[b],,joins[[a]],],c(1,3),sum),1,mean)
    lwr[,a]<-apply(apply(x.stored[data.years[b],,joins[[a]],],c(1,3),sum),1,quantile,probs=lwrq,na.rm=T)
    upr[,a]<-apply(apply(x.stored[data.years[b],,joins[[a]],],c(1,3),sum),1,quantile,probs=uprq,na.rm=T)
  }
  d.mid<-matrix(NA,ags,length(joins))
  d.upr<-matrix(NA,ags,length(joins))
  d.lwr<-matrix(NA,ags,length(joins))
  for (a in 1:length(joins)) {
    pos<-apply(d[data.years[b],1:ags,joins[[a]]],1,sum)
    tot<-apply(d[data.years[b],1:ags,],1,sum)
    for (k in 1:ags) {
      if (!is.nan(tot[k]) && tot[k]>0) {
        cis<-binom.test(pos[k],tot[k],conf.level=1-2*lwrq)$conf.int
        d.mid[k,a]<-pos[k]/tot[k]
        d.lwr[k,a]<-cis[1]
        d.upr[k,a]<-cis[2]
      }
    }
  }
  for (a in 1:length(joins)) {
    plot(1:ags,t="n",xaxt="n",xlim=c(1,ags),ylim=c(0,1),main=paste("Model fit for",disease.status[a],"in",2002+data.years[b]),xlab="Birth year",ylab="Prevalence")
    axis(1,1:ags,birth.years)
    wh<-which(!is.na(lwr[,a]))
    polygon(c(wh,rev(wh)),c(lwr[wh,a],rev(upr[wh,a])),col=a.col[a],border=NA)
    lines(1:ags,mid[,a],col=cols[a],lwd=2)
    segments(1:ags,d.lwr[,a],1:ags,d.upr[,a],col=a.col[a])
    points(1:ags,d.mid[,a],col=cols[a],pch=4)
  }
}
#
# Make predictions
#
fut<-3
x<-prediction(yrs=fut)
d<-add.data()
#
#  Plot predictions of the 3 diseases
#
mid<-matrix(NA,ags-fut,length(joins))
lwr<-matrix(NA,ags-fut,length(joins))
upr<-matrix(NA,ags-fut,length(joins))
for (a in 1:length(joins)) {
  mid[,a]<-apply(apply(x[tps,,joins[[a]],],c(1,3),sum),1,mean)
  lwr[,a]<-apply(apply(x[tps,,joins[[a]],],c(1,3),sum),1,quantile,probs=lwrq)
  upr[,a]<-apply(apply(x[tps,,joins[[a]],],c(1,3),sum),1,quantile,probs=uprq)
}
d.mid<-matrix(NA,ags-fut,length(joins))
d.upr<-matrix(NA,ags-fut,length(joins))
d.lwr<-matrix(NA,ags-fut,length(joins))
for (a in 1:length(joins)) {
  pos<-apply(d[tps,1:(ags-fut),joins[[a]]],1,sum)
  tot<-apply(d[tps,1:(ags-fut),],1,sum)
  for (k in 1:(ags-fut)) {
    if (!is.nan(tot[k]) && tot[k]>0) {
      cis<-binom.test(pos[k],tot[k],conf.level=1-2*lwrq)$conf.int
      d.mid[k,a]<-pos[k]/tot[k]
      d.lwr[k,a]<-cis[1]
      d.upr[k,a]<-cis[2]
    }
  }
}
#
# Make prediction plots
#
for (a in 1:length(joins)) {
  plot(1:(ags-fut),t="n",xaxt="n",xlim=c(1,ags-fut),ylim=c(0,1),main=paste("Predictions for",disease.status[a],"in 2016"),xlab="Birth year",ylab="Prevalence")
  axis(1,1:(ags-fut),birth.years[1:(ags-fut)])
  polygon(c(1:(ags-fut),rev(1:(ags-fut))),c(lwr[,a],rev(upr[,a])),col=a.col[a],border=NA)
  lines(1:(ags-fut),mid[,a],col=cols[a],lwd=2)
  segments(1:(ags-fut),d.lwr[1:(ags-fut),a],1:(ags-fut),d.upr[1:(ags-fut),a],col=a.col[a])
  points(1:(ags-fut),d.mid[1:(ags-fut),a],col=cols[a],pch=4)
}
dev.off()