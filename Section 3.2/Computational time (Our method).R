##### Premelianry, including source and library
#Enable command line arguments
#args = commandArgs(TRUE)
#p = as.numeric(args[1])
#amplitude_value = as.numeric(args[2])
#p = 50
#print(p)

source('KnockoffScreen_AL_vMar16.R')
source('ESGWAS.R')
source('GhostKnockoff_hello.R')

#source('GhostKnockoff.R')
#source('KnockoffScreen_AL_vMar16.R')
#Load SKAT haplotypes
#library(GhostKnockoff)
library(SKAT)
library(Matrix)
library(numDeriv)
library(glmnet)
library(corpcor)
#library(lassosum)
library(susieR)
library(CVXR)
#library(ghostbasil)
#library(Rfast)
###### Lasso
library(glmnet)

for(p in c(50,100,200,500)){
  ##### Set experimental setting, including n, p, k, and Sigma
  
  n = 1000         # number of observations
  #p = 200        # number of variables
  k = 5           # number of variables with nonzero coefficients
  
  amplitude <- 7
  rho_factor<-1
  
  
  
  ##### Experiment
  
  tau_1<-function(x)
  {
    l<-sort(x,decreasing = T)
    return(l[1]-median(l[-1]))
  }
  
  
  
  
  
  print(p)
  mu = rep(0,p)
  rho = 0.25*rho_factor
  Sigma = toeplitz(rho^(0:(p-1)))
  diag(Sigma)=1
  Sigma_chol<-chol(Sigma)
  
  
  
  
  
  FWER_1_d<-NULL
  PoWER_1_d<-NULL
  Time_spent<-NULL
  
  
  for (i in 1:200)
  { 
    time_0<-Sys.time()
    M<-19
    ####### Data generation
    mu = rep(0,p)
    rho = 0.25*rho_factor
    Sigma = toeplitz(rho^(0:(p-1)))
    diag(Sigma)=1
    beta_signal<-rep(0,p)
    Signal<-sample(1:p,k)
    beta_signal[Signal]<-sample(c(-1,1),k,replace = T)*(amplitude)/sqrt(n)
    # beta_signal[Signal]<-(amplitude)
    
    X<-matrix(rnorm(n*p),n,p)%*%Sigma_chol
    Y<-(X%*%beta_signal)[,1]+rnorm(n)
    X_mean<-colMeans(X)
    X<-t(X)-X_mean
    X_sd<-sqrt(rowMeans(X^2))
    X<-t(X/X_sd)
    Y<-(Y-mean(Y))
    Y<-Y/sqrt(mean(Y^2))
    T_0<-colSums(X*Y)/sqrt(n)
    time_1<-Sys.time()
    
    #### Optimization
    s<-create.solve_sdp_M(Sigma,M=M)
    diag_s<-Matrix(diag(s,length(s)))
    time_2<-Sys.time()
    
    #### Precalculation
    
    SigmaInv<-solve(Sigma)
    PP<-(diag(1,p)-diag_s%*%SigmaInv)
    time_3<-Sys.time()
    
    V1<-(M+1)/M*diag_s-diag_s%*%SigmaInv%*%diag_s
    V1<-(V1+t(V1))/2
    V1.right<-chol(V1)
    V2.right<-chol(diag_s)
    time_4<-Sys.time()
    
    M=19
    e1_tilde<-(t(V1.right)%*%rnorm(p))%x%t(rep(1,M))
    
    e2<-t(V2.right)%*%matrix(rnorm(p*M),p,M)
    
    e2_tilde<-e2-rowMeans(e2)
    
    time_5<-Sys.time()
    
    PZ<-(PP%*%T_0)%x%t(rep(1,M))
    T_null<-PZ+e1_tilde+e2_tilde
    
    time_6<-Sys.time()
    
    T_all_old<-cbind(T_0^2,T_null^2)
    kappa_old<-rowSums(T_null^2>T_0^2)
    tau_old<-apply(T_all_old,1,tau_1)
    ord_old<-order(tau_old,decreasing = T)
    kappa_old_ord<-kappa_old[ord_old]
    cut_old<-min(which(kappa_old_ord!=0))
    selected_old<-rep(0,p)
    selected_old[ord_old[(1:p)<cut_old]]<-1
    
    FWER_1_d<-c(FWER_1_d,sum(selected_old*(beta_signal==0))>=1)
    PoWER_1_d<-c(PoWER_1_d,sum(selected_old*(beta_signal!=0))/sum((beta_signal!=0)))
    
    
    
    time_7<-Sys.time()
    
    
    
    
    Time_one<-c(time_0,time_1,time_2,time_3,time_4,time_5,time_6,time_7)
    Time_spent<-rbind(Time_spent,Time_one[-1]-Time_one[-length(Time_one)])
    print(i)
  }
  ## Input the name of different columns
  colnames(Time_spent)<-c("data_generation","optimization","computing I-DSigma","cholesky decomposition","sampling","Computing Z","inference")
  ## Output 
  
  write.csv(Time_spent,paste0('time_efficient',p,'.csv'),row.names=F,col.names=T,quote=F,sep='\t')
}