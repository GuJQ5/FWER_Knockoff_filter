
source('KnockoffScreen_AL_vMar16.R')
source('ESGWAS.R')
source('GhostKnockoff_hello.R')
#SNR
p=50
pi <- rep(0,p)
n = 1000         # number of observations
#p = 200        # number of variables
k = 5           # number of variables with nonzero coefficients
#amplitude <- 7
rho_factor<-1
M<-1

#amplitude = 6 # signal amplitude (for noise level = 1)
mu = rep(0,p)
rho = 0.25*rho_factor
Sigma = toeplitz(rho^(0:(p-1)))
diag(Sigma)=1
Sigma_chol<-chol(Sigma)
amplitude<-3
##### Experiment
Time_spent<-NULL
FWER_old_d = NULL
PoWER_old_d<-NULL
M_deran = 50
eta = 0.99
for (i in 1:200)
{ 
  time_0<-Sys.time()
  M.d=1
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
  PP.d<-(diag(1,p)-diag_s%*%SigmaInv)
  time_3<-Sys.time()
  
  V1<-(M+1)/M*diag_s-diag_s%*%SigmaInv%*%diag_s
  V1<-(V1+t(V1))/2
  V1.right.d<-chol(V1)
  V2.right.d<-chol(diag_s)
  time_4<-Sys.time()
  
  pi <- rep(0,p)
  Time_iner<-NULL
  for (iter in 1:M_deran){ 
    time_4.0<-Sys.time()
    T_null_a<-(PP.d%*%T_0)%x%t(rep(1,M.d))
    T_null_b<-(t(V1.right.d)%*%rnorm(p))%x%t(rep(1,M.d))
    T_null_c<-t(V2.right.d)%*%matrix(rnorm(p*M.d),p,M.d)
    T_null_c<-T_null_c-rowMeans(T_null_c)
    T_null<-T_null_a+T_null_b+T_null_c
    ######################
    time_4.1<-Sys.time()
    T_all_old<-cbind(T_0^2,T_null^2)
    kappa_old<-rowSums(T_null^2>T_0^2)
    tau_old<-apply(T_all_old,1,tau_1)

    ord_old<-order(tau_old,decreasing = T)
    kappa_old_ord<-kappa_old[ord_old]
    cut_old<-min(which(kappa_old_ord!=0))
    
    pi[ord_old[(1:p)<cut_old]]<-pi[ord_old[(1:p)<cut_old]]+1
    time_4.2<-Sys.time()
    
    Time_iner1<-c(time_4.0,time_4.1,time_4.2)
    Time_iner<-rbind(Time_iner,Time_iner1[-1]-Time_iner1[-length(Time_iner1)])
    
  }
  
  time_5<-Sys.time()
  pi<-pi/M_deran
  selected_old_deran <- (pi>=eta)
  
  time_6<-Sys.time()
  
 
  
  
  Time_one<-c(time_0,time_1,time_2,time_3,time_4,time_5,time_6)
  Time_spent<-rbind(Time_spent,c(Time_one[-1]-Time_one[-length(Time_one)],colSums(Time_iner)))
 
  
  print(i)
  FWER_old_d<-c(FWER_old_d,sum(selected_old_deran*(beta_signal==0))>=1)
  PoWER_old_d<-c(PoWER_old_d,sum(selected_old_deran*(beta_signal!=0))/sum((beta_signal!=0)))
  

}


## Input the name of different columns
colnames(Time_spent)<-c("data_generation","optimization","computing I-DSigma","cholesky decomposition","loop inside","filnal filter","sampling","inference inside")
write.csv(Time_spent,paste0('time_deran',p,'.csv'),row.names=F,col.names=T,quote=F,sep='\t')


