#
# This file provides supporting functions for the MCMC.
# A key function is Dyson.full, which simulates the ode model.
#
library(mc2d)
library(mvtnorm)
library(far)
library(deSolve)
library(DescTools)
M<-orthonormalization(rep(1,8))
M<-t(M[,-1])
M1<-solve(rbind(M,rep(1,8)))
rwprior2<-function(v,k) {
  L<-length(v)
  return(k*sum((v[1:(L-2)]-2*v[2:(L-1)]+v[3:L])^2)/2)
}
rwprior1<-function(v,k) {
  L<-length(v)
  return(k*sum((v[1:(L-1)]-v[2:L])^2)/2)
}
projection<-function(x,tps,lambda,delta,mu,sigma) {
  y<-rep(x,tps)
  for (t in 2:tps) {
    y[t]<-(y[t-1]*(1-delta[t-1]-mu[t-1]-sigma[t-1])+(1-y[t-1])*(1-exp(-lambda[t-1])))/(1-mu[t-1]-delta[t-1]*y[t-1])
  }
  return(y)
}
adapt<-function(cur,accept,reject,target=0.45,a.min=0.15,a.max=0.5) {
  arate<-accept/(accept+reject)
  if ((arate<a.min || arate>a.max) && accept+reject>=adapt.freq) {
    if (accept==0) {
      sigma<-cur/2
    } else if (reject==0) {
      sigma<-cur*2
    } else {
      sigma<-cur*qnorm(target/2)/qnorm(arate/2)
    }
    return(c(sigma,0,0))
  } else {
    return(c(cur,accept,reject))
  }
}
# hcv hiv hsv
#1:000
#2:001
#3:010
#4:011
#5:100
#6:101
#7:110
#8:111
dyson.full<-function(t,p,parameters) {
  with(as.list(c(p,parameters)),{
    dp<-rep(NA,8)  
    hcv<-sum(p[5:8])
    hiv<-sum(p[c(3,4,7,8)])
    dp[1]<-p[1]*(delta[1]*hcv+delta[2]*hiv-lambda[1]-lambda[2]-lambda[3])
    dp[2]<-p[2]*(delta[1]*hcv+delta[2]*hiv-t[5]*lambda[1]-t[6]*lambda[2])          + p[1]*lambda[3]
    dp[3]<-p[3]*(delta[1]*hcv+delta[2]*hiv-t[3]*lambda[1]-t[4]*lambda[3]-delta[2]) + p[1]*lambda[2]
    dp[4]<-p[4]*(delta[1]*hcv+delta[2]*hiv-t[9]*lambda[1]-delta[2])                + p[2]*t[6]*lambda[2]+p[3]*t[4]*lambda[3]
    dp[5]<-p[5]*(delta[1]*hcv+delta[2]*hiv-t[1]*lambda[2]-t[2]*lambda[3]-delta[1]) + p[1]*lambda[1]
    dp[6]<-p[6]*(delta[1]*hcv+delta[2]*hiv-t[8]*lambda[2]-delta[1])                + p[2]*t[5]*lambda[1]+p[5]*t[2]*lambda[3]
    dp[7]<-p[7]*(delta[1]*hcv+delta[2]*hiv-t[7]*lambda[3]-delta[1]-delta[2])       + p[3]*t[3]*lambda[1]+p[5]*t[1]*lambda[2]
    dp[8]<-p[8]*(delta[1]*hcv+delta[2]*hiv-delta[1]-delta[2])                      + p[4]*t[9]*lambda[1]+p[6]*t[8]*lambda[2]+p[7]*t[7]*lambda[3]
    list(dp)
  })
}
plot.ode<-function(out) {
  plot(out[,1],t="n",ylim=0:1,xlab="Time",ylab="Proportion",xlim=c(0,max(out[,1]))) 
  for (i in 1:8) {
    lines(out[,1],apply(matrix(out[,2:(i+1)],length(out[,1]),i),1,sum))
  }
  lines(out[,1],apply(out[,1+5:8],1,sum),col="blue")
  lines(out[,1],apply(out[,1+c(3,4,7,8)],1,sum),col="red")
  lines(out[,1],apply(out[,1+2*1:4],1,sum),col="purple")
}
yearbyyear.full<-function(p,parameters,samplesPerYear=1) {
  tps<-length(parameters$lambda)/3 # NB this is usually tps-1
  times<-seq(0,1,by=1/samplesPerYear)
  out<-matrix(c(0,p),1,9)
  for (j in 1:tps) {
    parms=list(t=parameters$t,delta=parameters$delta[,j],lambda=parameters$lambda[,j])
    if (j==1) {
      y<-p
    } else {
      y<-out[1+samplesPerYear*(j-1),2:9]
    }
    temp<-ode(y=y, times=times,func=dyson.full, parms=parms,method="rk4")
    temp[,1]<-temp[,1]+j-1
    temp<-temp[-1,]
    out<-rbind(out,temp)
  }
  return(out)
}
load.data<-function(gender=1,race=1,byr.lwr=c(-Inf,1951,1963,1977,Inf)) {# 0-1950,1951-1962,1963-1976,1977-present
  if (gender=="men" | gender=="male") {gender<-1}
  if (gender=="women" | gender=="female") {gender<-2}
  if (race=="black") {race<-1}
  if (race=="white") {race<-2}
  dat<-read.csv("dataRANDOMIZED.csv")
  dat$byr<-dat$year-dat$age
  dat<<-dat
  gender<<-gender
  race<<-race
  tps<-11
  yrs<-2002+1:tps
  ags<-length(byr.lwr)-1
  # delete records with NA
  #wh<-which(!is.na(dat$hcv) & !is.na(dat$hsv))
  #dat<-dat[wh,]
  bin<-cbind(c(0,0,0,0,1,1,1,1),c(0,0,1,1,0,0,1,1),c(0,1,0,1,0,1,0,1))
  d<-array(NA,c(tps,ags,8))
  delta<-array(NA,c(3,tps-1,ags-1))
  delta[1,,]<-0.05868
  for (j in 1:tps) {
    if (j<tps) {delta[1,j,which(yrs[j]-byr.lwr[1:(ags-1)]<38)]<-0.016}
    for (k in 1:ags) {
      for (i in 1:8) {
        wh<-which(dat$gender==gender & dat$race==race & dat$hcv==bin[i,1] & dat$hiv==bin[i,2] & dat$hsv==bin[i,3] & dat$byr>=byr.lwr[k] & dat$byr<byr.lwr[k+1] & dat$year==yrs[j])
        d[j,k,i]<-length(wh)
      }
      if (j<tps & k<ags) {delta[2,j,k]<-get.delta(yrs[j],byr.lwr[k:(k+1)])}
    }
  }
  delta<<-delta
  return(d)
}
add.data<-function() {
  dat<-read.csv("data2016RANDOMIZED.csv")
  tps<<-tps+3
  ags<<-ags+3
  dat$byr<-2016-dat$age
  byrs<<-c(byrs[1:(length(byrs)-1)],1996,1997,1998,Inf)
  y18<<-pmax(1,byrs-1998+tps)[1:ags]
  yrs<-2002+1:tps
  d.old<<-d
  delta.old<-delta
  delta<-array(NA,c(3,tps-1,ags-1))
  delta[1,,]<-0.05868
  for (j in 1:(tps-1)) {
    delta[1,j,which(yrs[j]-byrs[1:(ags-1)]<38)]<-0.016
  }
  delta[2,1:(tps-4),1:(ags-4)]<-delta.old[2,,]
  d<-array(0,c(tps,ags,8))
  d[1:(tps-3),1:(ags-3),]<-d.old
  bin<-cbind(c(0,0,0,0,1,1,1,1),c(0,0,1,1,0,0,1,1),c(0,1,0,1,0,1,0,1))
  for (k in 1:ags) {
    for (i in 1:8) {
      wh<-which(dat$gender==gender & dat$race==race & dat$hcv==bin[i,1] & dat$hiv==bin[i,2] & dat$hsv==bin[i,3] & dat$byr>=byrs[k] & dat$byr<byrs[k+1])
      d[tps,k,i]<-length(wh)
    }
    if (k<ags) {
      for (j in 1:(tps-1)) {
        delta[2,j,k]<-get.delta(yrs[j],byrs[k:(k+1)])
      }
    }
  }
  delta<<-delta 
  return(d)
}
get.delta<-function(yr,byr.lwr) {
  if (yr<2006) {
    if (byr.lwr[1]==-Inf) {
      delta<-0.0308
    } else {
      if (byr.lwr[2]>yr-18) {byr.lwr[2]<-yr-18}
      yrs<-(byr.lwr[1]+1):byr.lwr[2]
      delta<-(length(intersect(yrs,1970:2015))*0.0107+length(intersect(yrs,1960:1969))*0.0197+length(intersect(yrs,1950:1959))*0.026+length(which(yrs<1950))*0.0308)/length(yrs)
    }          
  } else {
    if (byr.lwr[1]==-Inf) {
      delta<-0.0274
    } else {
      if (byr.lwr[2]>yr-18) {byr.lwr[2]<-yr-18}
      yrs<-(byr.lwr[1]+1):byr.lwr[2]
      delta<-(length(intersect(yrs,1973:2015))*0.011+length(intersect(yrs,1963:1972))*0.0125+length(intersect(yrs,1953:1962))*0.0181+length(which(yrs<1953))*0.0274)/length(yrs)
    }
  }
  return(delta)
}
x.plot<-function(k,main="",joins=list(5:8,c(3,4,7,8),2*1:4),a.col=c("#00FF0066","#FF000066","#0000FF66"),col=c("green","red","blue"),eps=0.025,conf.level=0.9, method="binomial") {
  sps<-dim(x.stored)[4]
  yrs<-2002+1:tps
  plot(yrs,t="n",ylim=c(0,1),xlim=c(2003,2002+tps),xlab="Year",ylab="Prevalence",main=main)
  j.d<-matrix(NA,tps,length(joins))
  j.d.not<-matrix(NA,tps,length(joins))
  mid<-matrix(NA,tps,length(joins))
  for (a in 1:length(joins)) {
    if (length(joins[[a]])==1) {
      j.x<-x.stored[,k,joins[[a]],]
      j.d[,a]<-d[,k,joins[[a]]]
    } else {
      j.x<-apply(x.stored[,k,joins[[a]],],c(1,3),sum)
      j.d[,a]<-apply(d[,k,joins[[a]]],1,sum)
    }
    j.d.not[,a]<-apply(d[,k,setdiff(1:8,joins[[a]])],1,sum)
    upr<-rep(NA,tps)
    lwr<-rep(NA,tps)
    for (j in y18[k]:tps) {
      mid[j,a]<-mean(j.x[j,])
      upr[j]<-quantile(j.x[j,],probs=uprq)
      lwr[j]<-quantile(j.x[j,],probs=lwrq)
    }
    polygon(c(yrs,rev(yrs)),c(upr,rev(lwr)),col=a.col[a],border=NA)
  }
  for (a in 1:length(joins)) {
    lines(yrs,mid[,a],col=col[a],lwd=2)
  }
  for (j in y18[k]:tps) {
    if (method=="binomial") { 
      if (!is.nan(sum(j.d[j,])) && j.d[j,1]+j.d.not[j,1]>0) {
        for (a in 1:length(joins)) {
          cis<-binom.test(j.d[j,a],j.d[j,a]+j.d.not[j,a],conf.level=conf.level)$conf.int
          segments(yrs[j]-(length(joins)/2+1/2)*eps+eps*a,cis[1],yrs[j]-(length(joins)/2+1/2)*eps+eps*a,cis[2],col=a.col[a])
          points(yrs[j]-(length(joins)/2+1/2)*eps+eps*a,j.d[j,a]/(j.d[j,a]+j.d.not[j,a]),col=col[a],pch=4)
        }
      }
    } else if (method=="goodman") { 
      if (!is.nan(sum(j.d[j,])) && sum(j.d[j,])>0) {
        cis<-MultinomCI(j.d[j,],conf.level=conf.level,method="goodman")
        for (a in 1:length(joins)) {
          segments(yrs[j]-(length(joins)/2+1/2)*eps+eps*a,cis[a,2],yrs[j]-(length(joins)/2+1/2)*eps+eps*a,cis[a,3],col=a.col[a])
          points(yrs[j]-(length(joins)/2+1/2)*eps+eps*a,cis[a,1],col=col[a],pch=4)
        }
      }
    } else {
      if (j.d[j,1]+j.d.not[j,1]>0) {
        for (a in 1:length(joins)) {
          points(yrs[j]-(length(joins)/2+1/2)*eps+eps*a,j.d[j,a]/(j.d[j,a]+j.d.not[j,a]),col=col[a],pch=4)
        }
      }
    }
  }
}
prediction<-function(yrs=3,usePrior=TRUE) {
  sps<-dim(x.stored)[4]
  x<-array(NA,c(tps+yrs,ags,8,sps))
  x[1:tps,,,]<-x.stored
  for (it in 1:sps) {
    for (k in 1:ags) {
      #params<-list(t=theta[1:ths],delta=matrix(delta[,y18[k]:(tps-1),k],3,tps-y18[k]),lambda=matrix(lambda[,y18[k]:(tps-1),k],3,tps-y18[k]))
      #x.proposed[y18[k]:tps,]<-yearbyyear(x[y18[k],k,],params)[,2:9]
      delta.mat<-matrix(delta[,tps-1,min(k,ags-1)],3,yrs)
      lambda.mat<-matrix(lambda.stored[,tps-1,min(k,ags-1),it],3,yrs) # copy params from earlier ages
      if (usePrior) { # draw params from prior
        if (k<ags) {
          lambda.mat[,1]<-rnorm(3,lambda.stored[,tps-1,k,it],k.time.stored[,it]^-0.5)
        } else { # add prior uncertainty for age; value borrowed from 1 age lower
          lambda.mat[,1]<-rnorm(3,lambda.stored[,tps-1,ags-1,it],k.age.stored[,it]^-0.5) 
          lambda.mat[,1]<-rnorm(3,lambda.mat[,1],k.time.stored[,it]^-0.5) # add prior uncertainty for time
        }
        for (j in 2:yrs) {
          lambda.mat[,j]<-rnorm(3,lambda.mat[,j-1],k.time.stored[,it]^-0.5) # add prior uncertainty for time
        }
      }
      params<-list(t=theta.stored[1:ths,it],delta=delta.mat,lambda=lambda.mat)
      x[tps:(tps+yrs),k,,it]<-yearbyyear(x.stored[tps,k,,it],params)[,2:9]
    }
  }
  return(x)
}

                                
