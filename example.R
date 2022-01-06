##################################
# Multimodal-integration example

# High-dimensional exposure and high-dimensional mediator
##################################

##################################
library("mvtnorm")

rm(list=ls())

# load functions
source("HD-CauseMediation.R")

# load model parameters
load("par.RData")
##################################

##################################
# Generate data

n<-100       # sample size n=100

r<-100       # # of exposure
p<-100       # # of mediator

set.seed(100)

# X
X<-rmvnorm(n,mean=rep(0,r),sigma=Sigma.X)
X.t<-X%*%Gamma[,1:q]

# M
e.M<-rmvnorm(n,mean=rep(0,p),sigma=Sigma.M)
M<-X.t%*%alpha+e.M

# Y
Y<-X.t%*%gamma+M%*%beta+rnorm(n,mean=0,sd=sigma.Y)
##################################

#################################
# method parameters

var.prop<-0.8          # PCA % of variance explained cut-off

standardize<-TRUE

max.itr<-10000
tol<-1e-6
trace<-FALSE

rho<-1
rho.increase<-FALSE

phi<-2                       # c0 tuning parameter

delta<-c(0.1,0.5,1,2)        # c1 tuning parameter
pi<-c(0.75,0.5,0.25,0)       # pi/(1-pi) = lambda1/lambda2    tuning the ratio

pi.delta.mat<-cbind(pi=rep(pi,each=length(delta)),delta=rep(delta,length(pi)))
pi.delta.mat<-pi.delta.mat[-which(pi.delta.mat[,"pi"]==0)[-1],]

lambda<-10^c(seq(-3,1,length.out=11),seq(1,2,length.out=3)[-1])          # lambda1 and lambda3
#################################

##################################
# run

# method: choose q based on variance proportion

re<-vector("list",nrow(pi.delta.mat))

for(ss in 1:nrow(pi.delta.mat))
{
  re[[ss]]<-vector("list",length(lambda))
  
  for(i in length(lambda):1)
  {
    if(i==length(lambda))
    {
      try(re[[ss]][[i]]<-HDCauseMediationPCA(X,M,Y,adaptive=TRUE,var.prop=var.prop,n.pc=q,lambda=lambda[i],pi=pi.delta.mat[ss,1],phi=phi,delta=pi.delta.mat[ss,2],rho=rho,
                                             standardize=standardize,max.itr=max.itr,tol=tol,rho.increase=rho.increase,trace=trace))
    }else
    {
      if(is.null(re[[ss]][[i+1]])==FALSE)
      {
        try(re[[ss]][[i]]<-HDCauseMediationPCA(X,M,Y,adaptive=TRUE,var.prop=var.prop,n.pc=q,lambda=lambda[i],pi=pi.delta.mat[ss,1],phi=phi,delta=pi.delta.mat[ss,2],rho=rho,
                                               standardize=standardize,max.itr=max.itr,tol=tol,rho.increase=rho.increase,trace=trace,
                                               alpha0=re[[ss]][[i+1]]$out.scaled$alpha,beta0=re[[ss]][[i+1]]$out.scaled$beta,gamma0=re[[ss]][[i+1]]$out.scaled$gamma))
      }else
      {
        try(re[[ss]][[i]]<-HDCauseMediationPCA(X,M,Y,adaptive=TRUE,var.prop=var.prop,n.pc=q,lambda=lambda[i],pi=pi.delta.mat[ss,1],phi=phi,delta=pi.delta.mat[ss,2],rho=rho,
                                               standardize=standardize,max.itr=max.itr,tol=tol,rho.increase=rho.increase,trace=trace))
      }
    }
    print(paste0("pi/delta index ",ss,", lambda index ",i," (",format((length(lambda)+1-i)/length(lambda)*100,digits=1,nsmall=1),"% complete)"))
  }
}

# BIC
zero.thred<-1e-3
BIC<-matrix(NA,length(lambda),nrow(pi.delta.mat))
for(ss in 1:nrow(pi.delta.mat))
{
  for(i in 1:length(lambda))
  {
    if(is.null(re[[ss]][[i]])==FALSE)
    {
      re.tmp<-re[[ss]][[i]]
      
      BIC[i,ss]<-log(n)*(sum(abs(re.tmp$out.scaled$IE)>zero.thred)+sum(abs(re.tmp$out.scaled$gamma)>zero.thred))-2*(-as.numeric(re.tmp$out.scaled$logLik[1]))
    }
  }
}
##################################

save.image("example.RData")




