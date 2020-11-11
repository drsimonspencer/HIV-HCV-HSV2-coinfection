#
# This function produces plots of the data (Figure 1).
# Four short runs of hiv.R (once for each gender/race combination) are needed to produce the required data files.
#
source('hivFunctions.R')
cols<-rev(c("#000000","#444444","#888888","#BBBBBB"))
pchs<-1:4
disease.status<-c("HCV","HIV","HSV2")
joins<-list(5:8,c(3,4,7,8),2*1:4)
data.years<-c(14,11,5,1)
group.names<-c("18-24","25-34","35-44","45-54","55+")
maxs<-c(0.5469613,0.3151515,0.8397626)
filenames<-c("Black females","White females","Black males","White males")   
#
pdf("fig_data.pdf",paper="a4",width=9.5,height=12,pointsize=10)
par(mfrow=c(4,3))
for (f in 1:4) {
  if (f==1) {
    filename<-"blackfemale"
  } else if (f==2) {
    filename<-"whitefemale"
  } else if (f==3) {
    filename<-"blackmale"
  } else {
    filename<-"whitemale"
  }
  load(paste0(filename,".Rdata"))
  d<-add.data()
  birth.years<-c(paste0("-",byrs[2]-1),byrs[2:ags])
  for (i in 1:3) {
    pos<-array(NA,c(3,4,5))
    tot<-array(NA,c(3,4,5))
    for (b in 1:length(data.years)) {
      for (a in 1:length(joins)) {
        if (b==1) {
          groups<-list(51:57,41:50,31:40,21:30,1:20)
        } else if (b==2) {
          groups<-list(48:54,38:47,28:37,18:27,1:17)
        } else if (b==3) {
          groups<-list(42:48,32:41,22:31,12:21,1:11)
        } else {
          groups<-list(38:44,28:37,18:27,8:17,1:7)
        }
        for (g in 1:length(groups)) { 
          pos[a,5-b,g]<-sum(d[data.years[b],groups[[g]],joins[[a]]])
          tot[a,5-b,g]<-sum(d[data.years[b],groups[[g]],])
        }
      }
    }
    p<-pos/tot
    plot(1:5,t="n",ylab="Prevalence",ylim=c(0,maxs[i]),main=disease.status[i],xaxt="n",xlab="")
    axis(1,1:5,group.names)
    if (i==2) {mtext(filenames[f],line=2.8,cex=1.1)}
    for (b in 1:length(data.years)) {
      lines(1:5,p[i,b,],col=cols[b])
      points(1:5,p[i,b,],pch=pchs[b],col=cols[b])
    }
    if (f==1 & i==1) {legend("topleft",legend=rev(data.years+2002),col=cols,pch=pchs,lty=1)}
    cat(f,i,max(p[i,,]),"\n")
  }
}
dev.off()

      