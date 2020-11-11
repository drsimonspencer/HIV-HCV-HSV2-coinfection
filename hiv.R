#
# This file runs the MCMC - shortened version for testing takes about 7 mins.
# Uncomment one line between 12-15 to run the different datasets. Comment line 45 for full length run (takes considerable time ie days).
# Doing a second run loads the covariance matrix from the previous run, providing greatly improved mixing.
source('hivFunctions.R')
yearbyyear<-yearbyyear.full
tps<-11 # timepoints
ags<-54; byrs<-c(-Inf,1943:1995,Inf) # ages
#ags<-14; byrs<-c(-Inf,1983:1995,Inf)
y18<-pmax(1,byrs-1995+tps)[1:ags]
ths<-9
d<-load.data(gender=1,byr.lwr=byrs);filename<-"blackmale" # load data
#d<-load.data(gender=2,byr.lwr=byrs);filename<-"blackfemale"
#d<-load.data(gender=1,race=2,byr.lwr=byrs);filename<-"whitemale"
#d<-load.data(gender=2,race=2,byr.lwr=byrs);filename<-"whitefemale"
#d<-add.data() uncomment this line to fit to all datasets
simulated<-FALSE
if (simulated) { # simulate data, overwriting real data
  lambda.true<-array(rep(0.01+10:1/200,each=3),c(3,tps-1,ags)) # true infection rate  (3xtps-1)xags
  delta<-array(rep(c(0.0345,0.020,NA),ags*(tps-1)),c(3,tps-1,ags)) # known death rate of infected individuals (minus natural death rate)
  theta.true<-c(2,1,1,1,1,5,1,1,1)
  size<-matrix(50,tps,ags)#c(50,0,0,0,0,50,0,0,0,0,50) # number of samples tested
  x.true<-array(rep(c(0.79,0.03,0.03,0.03,0.03,0.03,0.03,0.03),each=tps*ags),c(tps,ags,8))
  d<-array(NA,c(tps,ags,8))
  for (k in 1:ags) {
    params<-list(t=theta.true,delta=delta[,,k],lambda=lambda.true[,,k])
    x.true[,k,]<-yearbyyear(x.true[1,k,],params)[,2:9]
    for (j in 1:tps) {
      d[j,k,]<-rmultinom(1,size=size[j,k],x.true[j,k,])
    }
  }
}
set.seed(102)
# priors
a.time<-1
b.time<-0.01
a.age<-1
b.age<-0.01
max.lambda<-Inf
alpha.init<-rep(0.8/8,8)
a.theta<-rep(1,ths)
b.theta<-rep(1,ths)
# initialise MCMC
iters<-120000;thinning<-50;burnin<-20000;adapt.freq<-500 # These are the values used in the paper.
iters<-60;thinning<-1;burnin<-10;adapt.freq<-10 # short run to test code
blocks<-TRUE
sing<-TRUE
lambda<-array(NA,c(3,tps-1,ags-1))
for (k in 1:(ags-1)) { # last age-group has no lambda
  lambda[,y18[k]:(tps-1),k]<-0.005+0.005*rbeta(3*(tps-y18[k]),2,2)
}
lambda.not.na<-which(!is.na(lambda[1,,]))
n.lambdas<-sum(!is.na(lambda[1,,]))
a.lambda.mean<-100/n.lambdas # mean incidence is 1%
lambda.stored<-array(NA,c(3,tps-1,ags-1,(iters-burnin)/thinning))
lambda.accept<-array(0,c(3,tps-1,ags-1,2))
lambda.reject<-array(0,c(3,tps-1,ags-1,2))
lambda.sigma<-array(0.02,c(3,tps-1,ags-1))
k.time<-rep(10,3)
k.time.stored<-matrix(NA,3,(iters-burnin)/thinning)
k.age<-rep(10,3)
k.age.stored<-matrix(NA,3,(iters-burnin)/thinning)
x<-array(NA,c(tps,ags,8))
for (k in 1:ags) {
  x[y18[k]:tps,k,]<-rep(rdirichlet(1,alpha.init+d[y18[k],k,]),each=tps-y18[k]+1)
}
x.not.na<-which(!is.na(x) & d!=0)
theta<-c(rep(1,ths),rep(NA,3))
theta.stored<-matrix(NA,length(theta),(iters-burnin)/thinning)
theta.accept<-rep(0,ths)
theta.reject<-rep(0,ths)
theta.sigma<-rep(0.5,ths)
theta.accept.block<-0
theta.reject.block<-0
theta.stored.all<-matrix(NA,length(theta),iters)
adapt.factor<-NA  # artificial reduction to try to improve acceptance.
if (file.exists(paste0("Sigma",filename,".Rdata"))) { # load previous posterior covariance matrix if available
  load(paste0("Sigma",filename,".Rdata"))  
} else {
  theta.Sigma<-diag(rep(0.0001,length(theta)))
}
for (k in 1:(ags-1)) { # last age-group has no lambda
  params<-list(t=theta[1:ths],delta=matrix(delta[,y18[k]:(tps-1),k],3,tps-y18[k]),lambda=matrix(lambda[,y18[k]:(tps-1),k],3,tps-y18[k]))
  x[y18[k]:tps,k,]<-yearbyyear(x[y18[k],k,],params)[,2:9]
}
x.stored<-array(NA,c(tps,ags,8,(iters-burnin)/thinning))
x.accept<-rep(0,ags)
x.reject<-rep(0,ags)
x.sigma<-rep(20,ags)
# MCMC loop
st<-Sys.time()
for (it in 1:iters) {
  adapt.accept<-1+0.2*(100/(100+it))
  adapt.reject<-adapt.accept^(0.44/(0.44-1))
  # update lambda with rw
  for (i in 1:3) {
    for (k in 1:(ags-1)) { # this loop is embarassing (two apart)
      for (j in y18[k]:(tps-1)) {
        proposal<-lambda
        if (it%%2==1 && j>y18[k] && j<tps-1 && k>1 && k<ags-tps+j) {# Conditional prior proposal
          cpp<-TRUE 
          cpp.mu<-(k.time[i]*(lambda[i,j-1,k]+lambda[i,j+1,k])/2+k.age[i]*(lambda[i,j,k-1]+lambda[i,j,k+1])/2)/(k.time[i]+k.age[i])
          cpp.sd<-(2*k.time[i]+2*k.age[i])^(-1/2)
          proposal[i,j,k]<-rnorm(1,cpp.mu,cpp.sd)
        } else { # MH RW
          cpp<-FALSE
          proposal[i,j,k]<-rnorm(1,lambda[i,j,k],lambda.sigma[i,j,k])
        }
        if (proposal[i,j,k]>0 && proposal[i,j,k]<max.lambda) {
          params<-list(t=theta[1:ths],delta=matrix(delta[,j:(tps-1),k],3,tps-j),lambda=matrix(proposal[,j:(tps-1),k],3,tps-j))
          x.proposed<-x[,k,]
          x.proposed[j:tps,]<-yearbyyear(x[j,k,],params)[,2:9]
          if (min(x.proposed[j:tps,])>=0 && max(x.proposed[j:tps,])<=1) {
            ap<--a.lambda.mean*(proposal[i,j,k]-lambda[i,j,k])#/(tps-1)/ags
            if (!cpp) {
              if (k<ags-1) {
                ap<-ap+rwprior1(proposal[i,max(y18[k],j-1):min(tps-1,j+1),k],-k.time[i])-rwprior1(lambda[i,max(y18[k],j-1):min(tps-1,j+1),k],-k.time[i])
              }
              ap<-ap+rwprior1(proposal[i,j,max(1,k-1):min(ags-tps+j,k+1)],-k.age[i])-rwprior1(lambda[i,j,max(1,k-1):min(ags-tps+j,k+1)],-k.age[i])
            } # q cancels with prior for cpp
            #ap<-ap+sum(dmultinomial(d[(j+1):tps,k,],size=size[(j+1):tps,k],x.proposed[(j+1):tps,],log=TRUE)-dmultinomial(d[(j+1):tps,k,],size=size[(j+1):tps,k],x[(j+1):tps,k,],log=TRUE)) # alternative coding of likelihood
            ap<-ap+sum(d[(j+1):tps,k,]*(log(x.proposed[(j+1):tps,])-log(x[(j+1):tps,k,])))
            u<-runif(1)     
            if (u<=exp(ap)) {
              lambda[i,j,k]<-proposal[i,j,k]
              x[,k,]<-x.proposed
              lambda.accept[i,j,k,1+cpp]<-lambda.accept[i,j,k,1+cpp]+1
              if (!cpp) {lambda.sigma[i,j,k]<-lambda.sigma[i,j,k]*adapt.accept}
            } else {
              lambda.reject[i,j,k,1+cpp]<-lambda.reject[i,j,k,1+cpp]+1
              if (!cpp) {lambda.sigma[i,j,k]<-lambda.sigma[i,j,k]*adapt.reject}
            }
          } else {
            lambda.reject[i,j,k,1+cpp]<-lambda.reject[i,j,k,1+cpp]+1
            if (!cpp) {lambda.sigma[i,j,k]<-lambda.sigma[i,j,k]*adapt.reject}
          } 
        } else {
          lambda.reject[i,j,k,1+cpp]<-lambda.reject[i,j,k,1+cpp]+1
          if (!cpp) {lambda.sigma[i,j,k]<-lambda.sigma[i,j,k]*adapt.reject}
        }
      }
    }
  }
  for (i in 1:3) {
    # reset lambda means in theta
    temp<-lambda[i,,]
    theta[ths+i]<-mean(temp[lambda.not.na])
    # update k.time
    rwp<-0
    num.diffs<-0
    for (k in 1:(ags-2)) {
      rwp<-rwp+rwprior1(lambda[i,y18[k]:(tps-1),k],1)
      num.diffs<-num.diffs+tps-y18[k]-1
    }
    k.time[i]<-rgamma(1,a.time+num.diffs/2,rate=b.time+rwp)
    # update k.age
    rwp<-0
    num.diffs<-0
    for (j in 1:(tps-1)) {
      rwp<-rwp+rwprior1(lambda[i,j,1:(ags-tps+j)],1)
      num.diffs<-num.diffs+ags-tps+j-1
    }
    k.age[i]<-rgamma(1,a.age+num.diffs/2,rate=b.age+rwp)
  }
  # update initial values in x with MH for k<ags
  for (k in 1:(ags-1)) {
    x.proposed<-x[,k,]
    #x.proposed[y18[k],]<-M1%*%c(rmvnorm(1,M%*%x[y18[k],k,],x.sigma[k]*diag(7)),1)
    x.proposed[y18[k],]<-rdirichlet(1,alpha.init+d[y18[k],k,]+x.sigma[k]*x[y18[k],k,]) # This proposal would be the posterior if it wasn't for the other timepoints
    if (min(x.proposed[y18[k],])>=0 && max(x.proposed[y18[k],])<=1 && k<ags) {
      params<-list(t=theta[1:ths],delta=matrix(delta[,y18[k]:(tps-1),k],3,tps-y18[k]),lambda=matrix(lambda[,y18[k]:(tps-1),k],3,tps-y18[k]))
      x.proposed[y18[k]:tps,]<-yearbyyear(x.proposed[y18[k],],params)[,2:9]
    }
    if (min(x.proposed[y18[k]:tps,])>=0 && max(x.proposed[y18[k]:tps,])<=1) {
      ap<-sum(d[y18[k]:tps,k,]*(log(x.proposed[y18[k]:tps,])-log(x[y18[k]:tps,k,])))
      ap<-ap+sum((alpha.init-1)*(log(x.proposed[y18[k],])-log(x[y18[k],k,])))
      ap<-ap+log(ddirichlet(x[y18[k],k,],alpha.init+d[y18[k],k,]+x.sigma[k]*x.proposed[y18[k],]))
      ap<-ap-log(ddirichlet(x.proposed[y18[k],],alpha.init+d[y18[k],k,]+x.sigma[k]*x[y18[k],k,]))
      if (k==54) {cat(ap,"\n")}
      u<-runif(1)     
      if (u<=exp(ap)) {
        x[,k,]<-x.proposed
        x.accept[k]<-x.accept[k]+1
        if (it<burnin) {x.sigma[k]<-max(0,x.sigma[k]-3)}
      } else {
        x.reject[k]<-x.reject[k]+1
        if (it<burnin) {x.sigma[k]<-x.sigma[k]+1}
      }
    } else {
      x.reject[k]<-x.reject[k]+1
      if (it<burnin) {x.sigma[k]<-x.sigma[k]+1}
    }
  }
  # update x with Gibbs for k=ags
  x[y18[ags],ags,]<-rdirichlet(1,alpha.init+d[y18[ags],ags,])
  # update theta  #swap to log normal proposals?
  if (blocks) {# block updates for theta, including the mean of the lambdas
    proposal<-rmvnorm(1,c(theta[1:ths],0,0,0),theta.Sigma)
    if (min(proposal[1:ths])>0 && min(lambda[1,,],na.rm=T)+proposal[ths+1]>0 && min(lambda[2,,],na.rm=T)+proposal[ths+2]>0 && min(lambda[3,,],na.rm=T)+proposal[ths+3]>0 && max(lambda[1,,],na.rm=T)+proposal[ths+1]<max.lambda && max(lambda[2,,],na.rm=T)+proposal[ths+2]<max.lambda && max(lambda[3,,],na.rm=T)+proposal[ths+3]<max.lambda) {
    #if (min(proposal[1:ths])>0 && min(apply(lambda,1,min,na.rm=T)+proposal[(ths+1):(ths+3)])>0 && max(apply(lambda,1,max,na.rm=T)+proposal[(ths+1):(ths+3)]<max.lambda)) { 
      x.proposed<-x
      for (k in 1:(ags-1)) {
        params<-list(t=proposal[1:ths],delta=matrix(delta[,y18[k]:(tps-1),k],3,tps-y18[k]),lambda=matrix(lambda[,y18[k]:(tps-1),k],3,tps-y18[k])+rep(proposal[(ths+1):(ths+3)],tps-y18[k]))
        x.proposed[y18[k]:tps,k,]<-yearbyyear(x[y18[k],k,],params)[,2:9]
      }
      if (max(x.proposed[x.not.na])<=1 && min(x.proposed[x.not.na])>=0) {
        ap<-sum(dgamma(proposal[1:ths],a.theta[1:ths],b.theta[1:ths],log=TRUE)-dgamma(theta[1:ths],a.theta[1:ths],b.theta[1:ths],log=TRUE))
        ap<-ap-a.lambda.mean*sum(proposal[(ths+1):(ths+3)])*n.lambdas # proposal 5:7 is the difference between the proposal and theta.
        ap<-ap+sum(d[x.not.na]*(log(x.proposed[x.not.na])-log(x[x.not.na])))
        u<-runif(1)
        if (u<=exp(ap)) {
          theta[1:ths]<-proposal[1:ths]
          theta[(ths+1):(ths+3)]<-theta[(ths+1):(ths+3)]+proposal[(ths+1):(ths+3)]
          lambda[1,,]<-lambda[1,,]+proposal[ths+1]
          lambda[2,,]<-lambda[2,,]+proposal[ths+2]
          lambda[3,,]<-lambda[3,,]+proposal[ths+3]
          x<-x.proposed
          theta.accept.block<-theta.accept.block+1
        } else {
          theta.reject.block<-theta.reject.block+1
        }
      } else {
        theta.reject.block<-theta.reject.block+1
      }
    } else {
      theta.reject.block<-theta.reject.block+1       
    }
  }
  if (sing) {# single site updates for theta
    for (i in 1:ths) {
      proposal<-theta
      proposal[i]<-rnorm(1,theta[i],theta.sigma[i])
      if (proposal[i]>0) {
        x.proposed<-x
        for (k in 1:(ags-1)) {
          params<-list(t=proposal[1:ths],delta=matrix(delta[,y18[k]:(tps-1),k],3,tps-y18[k]),lambda=matrix(lambda[,y18[k]:(tps-1),k],3,tps-y18[k]))
          x.proposed[y18[k]:tps,k,]<-yearbyyear(x[y18[k],k,],params)[,2:9]
        }      
        if (max(x.proposed[x.not.na])<=1 && min(x.proposed[x.not.na])>=0) {
          ap<-dgamma(proposal[i],a.theta[i],b.theta[i],log=TRUE)-dgamma(theta[i],a.theta[i],b.theta[i],log=TRUE)
          ap<-ap+sum(d[x.not.na]*(log(x.proposed[x.not.na])-log(x[x.not.na])))
          u<-runif(1)
          if (u<=exp(ap)) {
            theta[i]<-proposal[i]
            x<-x.proposed
            theta.accept[i]<-theta.accept[i]+1
            theta.sigma[i]<-theta.sigma[i]*adapt.accept
          } else {
            theta.reject[i]<-theta.reject[i]+1
            theta.sigma[i]<-theta.sigma[i]*adapt.reject
          }
        } else {
          theta.reject[i]<-theta.reject[i]+1
          theta.sigma[i]<-theta.sigma[i]*adapt.reject
        }
      } else {
        theta.reject[i]<-theta.reject[i]+1
        theta.sigma[i]<-theta.sigma[i]*adapt.reject
      }
    }
  }      
  # take a sample
  if (it>burnin && it%%thinning==0) {
    s<-(it-burnin)%/%thinning
    lambda.stored[,,,s]<-lambda
    k.time.stored[,s]<-k.time
    k.age.stored[,s]<-k.age
    x.stored[,,,s]<-x
    theta.stored[,s]<-theta
  }
  theta.stored.all[,it]<-theta
  # adapt
  if (it%%adapt.freq==0) {
    cat(it,"iterations completed in",round(difftime(Sys.time(),st,units="mins"),2),"minutes. Block acceptance rate",
      round(theta.accept.block/(theta.accept.block+theta.reject.block),3),"adapt factor",adapt.factor,"\n")
    cat("theta.accept:",round(theta.accept/(theta.reject+theta.accept),3),"\n")
    cat("theta.sigma:",round(theta.sigma,3),"\n")
    #cat("sigma:",round(lambda.accept/(lambda.reject+lambda.accept),3),"\n")
    #cat("x:",round(x.accept/(x.reject+x.accept),3),"\n")
    # adapt lambda                                                                                                               
    #for (i in 1:3) {
    #  for (k in 1:(ags-1)) {
    #    for (j in y18[k]:(tps-1)) {
    #      temp<-adapt(lambda.sigma[i,j,k],lambda.accept[i,j,k,2],lambda.reject[i,j,k,2])
    #      lambda.sigma[i,j,k]<-temp[1]
    #      lambda.accept[i,j,k,2]<-temp[2]
    #      lambda.reject[i,j,k,2]<-temp[3]
    #      if (temp[2]==0 & temp[3]==0 & it>burnin) {cat("lambda",i,j,k,temp[1],"\n")}
    #    }
    #  }
    #}
    # adapt x is done earlier
    k<-which.min(x.accept[1:(ags-1)])[1]
    cat("Worst x ",k,x.sigma[k],round(x.accept[k]/(x.reject[k]+x.accept[k]),3),"\n")
    k<-which.max(x.accept[1:(ags-1)])[1]
    cat("Best x",k,x.sigma[k],round(x.accept[k]/(x.reject[k]+x.accept[k]),3),"\n")
    # adapt theta
    if (blocks) {
    #  temp<-adapt(sqrt(adapt.factor),theta.accept.block,theta.reject.block,0.1,0.05,0.4)
    #  adapt.factor<-temp[1]^2
    #  theta.accept.block<-temp[2]
    #  theta.reject.block<-temp[3]
    #  theta.Sigma<-adapt.factor*2.38^2/length(theta)*cov(t(theta.stored.all[,1:it]))#+0.01*diag(7)# artificial reduction to try to improve acceptance.
      theta.Sigma<-2.38^2/length(theta)*cov(t(theta.stored.all[,floor(it/2):it])) # probably best to stick with this.
    }
    #if (sing) {
    #  for (i in 1:ths) {
    #    temp<-adapt(theta.sigma[i],theta.accept[i],theta.reject[i])
    #    theta.sigma[i]<-temp[1]
    #    theta.accept[i]<-temp[2]
    #    theta.reject[i]<-temp[3]
    #    if (temp[2]==0 & temp[3]==0 & it>burnin) {cat("theta",i,temp[1],"\n")}
    #  }
    #}
  }
}
save.image(paste0(filename,".Rdata"))
save(theta.Sigma,file=paste0("Sigma",filename,".Rdata"))    
    
