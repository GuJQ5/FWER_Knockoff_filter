#source("GhostKnockoff.R")

##### Premelianry, including source and library

#Enable command line arguments


#amplitude_value = as.numeric(args[2])

args = commandArgs(TRUE)
 p = as.numeric(args[1])
 rho = as.numeric(args[2])

#p =500

#rho =0.5


 source('KnockoffScreen_AL_vMar16.R')
 source('ESGWAS.R')
 source('GhostKnockoff_hello.R')



#source('GhostKnockoff.R')
#source('KnockoffScreen_AL_vMar16.R')
#Load SKAT haplotypes
#library(GhostKnockoff)
library(SKAT)
library(Matrix)
library(MASS)
library(mvtnorm)
library(numDeriv)
library(glmnet)
library(corpcor)
#library(lassosum)
library(susieR)
#library(CVXR)
#library(ghostbasil)
#library

###### Lasso
#library(glmnet)






M_deran = 50
eta = 0.99
M_cderan = 10
eta_c = 0.89
v_list <- seq(1,20,by=1)
p_list <- pnbinom(1,v_list,.5,lower.tail = FALSE)+dnbinom(1,v_list,.5) 
if(sum(p_list<=0.05)>0){
  v0 <- max(which(p_list<=0.05))
  prob <- (0.05-p_list[v0])/(p_list[v0+1] - p_list[v0])
}else{
  v0 <- 0
  prob <- (0.05)/(p_list[1])  
}


tau_cal<-function(x)
{
  l<-sort(x,decreasing = T)
  return(l[1]-median(l[-1]))
}


tau_new_function<-function(x)
{
  l = abs(x-median(x))
  ls = sort(l,decreasing=T)
  
  return(ls[1])
}


n=500
#p=100
k=5
#rho=0.5
amp.seq = c(9)




FWER_old_deran<-NULL
PoWER_old_deran<-NULL


FWER_old_deran_sd<-NULL
PoWER_old_deran_sd<-NULL

FWER_old<-NULL
PoWER_old<-NULL

FWER_old_sd<-NULL
PoWER_old_sd<-NULL

FWER_new<-NULL
PoWER_new<-NULL

FWER_deran<-NULL
PoWER_deran<-NULL

FWER_deran_sd<-NULL
PoWER_deran_sd<-NULL

FWER_vanilla<-NULL
PoWER_vanilla<-NULL

FWER_vanilla_sd<-NULL
PoWER_vanilla_sd<-NULL


#Sigma<-toeplitz(rho^(1:p-1))
#diag(Sigma)<-1

###### Compound Symmetry
Sigma<-diag(1-rho,p)+rho


Sigma_chol<-chol(Sigma)
SigmaInv<-solve(Sigma)


M.d<-1
s<-create.solve_sdp_M(Sigma,M=M.d)
D<-Matrix(diag(s,length(s)))
D = 0.9*D
PP.d<-(diag(1,p)-D%*%SigmaInv)
V1.d<-(M.d+1)/M.d*D-D%*%SigmaInv%*%D
V1.d<-(V1.d+t(V1.d))/2
V1.right.d<-chol(V1.d)
V2.right.d<-chol(D)


M<-19
s<-create.solve_sdp_M(Sigma,M=M)
D<-Matrix(diag(s,length(s)))

PP<-(diag(1,p)-D%*%SigmaInv)
V1<-(M+1)/M*D-D%*%SigmaInv%*%D
V1<-(V1+t(V1))/2
V1.right<-chol(V1)
V2.right<-chol(D)



M.c<-5
s<-create.solve_sdp_M(Sigma,M=M.c)
D<-Matrix(diag(s,length(s)))

PP.c<-(diag(1,p)-D%*%SigmaInv)
V1.c<-(M.c+1)/M.c*D-D%*%SigmaInv%*%D
V1.c<-(V1.c+t(V1.c))/2
V1.right.c<-chol(V1.c)
V2.right.c<-chol(D)




for (amplitude in amp.seq){
  
  
  
  FWER_old_d<-NULL
  PoWER_old_d<-NULL
  
  
  FWER_old_a<-NULL
  PoWER_old_a<-NULL
  
  
  FWER_v<-NULL
  PoWER_v<-NULL
  
    
  FWER_new_d<-NULL
  PoWER_new_d<-NULL
  
 

  
  
  result.table.iner = NULL
  
  for(iter in 1:2000)
  {
    ####### Data generation
    beta_signal<-rep(0,p)
    Signal<-sample(1:p,k)
    beta_signal[Signal]<-sample(c(-1,1),k,replace = T)*(amplitude)/sqrt(n)
    # beta_signal[Signal]<-(amplitude)
    
    X<-matrix(rnorm(n*p),n,p)%*%Sigma_chol
    Y<-(X%*%beta_signal)[,1]+rnorm(n)
    ######################
    X_mean<-colMeans(X)
    X<-t(X)-X_mean
    X_sd<-sqrt(rowMeans(X^2))
    X<-t(X/X_sd)
    Y<-(Y-mean(Y))
    Y<-Y/sqrt(mean(Y^2))
    T_0<-colSums(X*Y)/sqrt(n)
    
    ##########################################derandomized##########################################
    ######################
    
    
    ######################
    
    
    
    M = 19

    v <- 1
    T_null_a<-(PP%*%T_0)%x%t(rep(1,M))
    T_null_b<-(t(V1.right)%*%rnorm(p))%x%t(rep(1,M))
    T_null_c<-t(V2.right)%*%matrix(rnorm(p*M),p,M)
    T_null_c<-T_null_c-rowMeans(T_null_c)
    T_null<-T_null_a+T_null_b+T_null_c
    T_all_old<-cbind(T_0^2,T_null^2)
    kappa_old<-rowSums(T_null^2>T_0^2)
    tau_old<-apply(T_all_old,1,tau_cal)
    ord_old<-order(tau_old,decreasing = T)
    kappa_old_ord<-kappa_old[ord_old]
    cut_old<-which(kappa_old_ord!=0)[v]
    cut_no<-which(kappa_old_ord!=0)[1:v-1]
    
    selected_old<-rep(0,p)
    selected_old[ord_old[(1:p)<cut_old]]<-1
    selected_old[ord_old[cut_no]]<-0
    
    FWER_old_a<-c(FWER_old_a,sum(selected_old*(beta_signal==0))>=1)
    PoWER_old_a<-c(PoWER_old_a,sum(selected_old*(beta_signal!=0))/sum((beta_signal!=0)))

    #####################
    ##################### 
    M.d = 1
    pi <- rep(0,p)
    
    for (m in 1:M_deran){
      
      T_null_a<-(PP.d%*%T_0)%x%t(rep(1,M.d))
      T_null_b<-(t(V1.right.d)%*%rnorm(p))%x%t(rep(1,M.d))
      T_null_c<-t(V2.right.d)%*%matrix(rnorm(p*M.d),p,M.d)
      T_null_c<-T_null_c-rowMeans(T_null_c)
      T_null<-T_null_a+T_null_b+T_null_c
      ######################
      T_all_old<-cbind(T_0^2,T_null^2)
      kappa_old<-rowSums(T_null^2>T_0^2)
      tau_old<-apply(T_all_old,1,tau_cal)
      #####################
      ord_old<-order(tau_old,decreasing = T)
      kappa_old_ord<-kappa_old[ord_old]
      cut_old<-min(which(kappa_old_ord!=0))
      
      pi[ord_old[(1:p)<cut_old]]<-pi[ord_old[(1:p)<cut_old]]+1
      #selected_old[ord_old[(1:p)<cut_old]]<-1
      #####################
      # ord_new<-order(tau_new,decreasing = T)
      # kappa_new_ord<-kappa_new[ord_new]
      # cut_new<-min(which(kappa_new_ord!=0))
      # #selected_new[ord_new[(1:p)<cut_new]]<-1
      # pi2[ord_new[(1:p)<cut_new]]=pi2[ord_new[(1:p)<cut_new]]+1
    }
    
    pi<-pi/M_deran
    selected_old_deran <- (pi>=eta)
    
    
    
    FWER_old_d<-c(FWER_old_d,sum(selected_old_deran*(beta_signal==0))>=1)
    PoWER_old_d<-c(PoWER_old_d,sum(selected_old_deran*(beta_signal!=0))/sum((beta_signal!=0)))
    
    
    
    
        #####################
    ##################### 
    M.c = 5
    pi.c <- rep(0,p)
    
    for (m in 1:M_cderan){
      
      T_null_a<-(PP.c%*%T_0)%x%t(rep(1,M.c))
      T_null_b<-(t(V1.right.c)%*%rnorm(p))%x%t(rep(1,M.c))
      T_null_c<-t(V2.right.c)%*%matrix(rnorm(p*M.c),p,M.c)
      T_null_c<-T_null_c-rowMeans(T_null_c)
      T_null<-T_null_a+T_null_b+T_null_c
      ######################
      T_all_old<-cbind(T_0^2,T_null^2)
      kappa_old<-rowSums(T_null^2>T_0^2)
      tau_old<-apply(T_all_old,1,tau_cal)
      #####################
      ord_old<-order(tau_old,decreasing = T)
      kappa_old_ord<-kappa_old[ord_old]
      cut_old<-min(which(kappa_old_ord!=0))
      
      pi[ord_old[(1:p)<cut_old]]<-pi[ord_old[(1:p)<cut_old]]+1
      #selected_old[ord_old[(1:p)<cut_old]]<-1
      #####################
      # ord_new<-order(tau_new,decreasing = T)
      # kappa_new_ord<-kappa_new[ord_new]
      # cut_new<-min(which(kappa_new_ord!=0))
      # #selected_new[ord_new[(1:p)<cut_new]]<-1
      # pi2[ord_new[(1:p)<cut_new]]=pi2[ord_new[(1:p)<cut_new]]+1
    }
    
    pi<-pi/M_cderan
    selected_old_deran <- (pi>=eta_c)
    
    
    
    FWER_new_d<-c(FWER_new_d,sum(selected_old_deran*(beta_signal==0))>=1)
    PoWER_new_d<-c(PoWER_new_d,sum(selected_old_deran*(beta_signal!=0))/sum((beta_signal!=0)))
    
    
    
    
    
    
    
    
    
    
    
    
    
    ## Run v-knockoffs 
    
 
   
    ######################
    M=1
    v_list <- seq(1,20,by=1)
    p_list <- pnbinom(1,v_list,1-1/(M+1),lower.tail = FALSE)+dnbinom(1,v_list,1-1/(M+1)) 
    if(sum(p_list<=0.05)>0){
      v0 <- max(which(p_list<=0.05))
      prob <- (0.05-p_list[v0])/(p_list[v0+1] - p_list[v0])
    }else{
      v0 <- 0
      prob <- (0.05)/(p_list[1])  
    }
    
    v <- floor(v0)+rbinom(1,1,prob)
    T_null_a<-(PP.d%*%T_0)%x%t(rep(1,M))
    T_null_b<-(t(V1.right.d)%*%rnorm(p))%x%t(rep(1,M.d))
    T_null_c<-t(V2.right.d)%*%matrix(rnorm(p*M),p,M.d)
    T_null_c<-T_null_c-rowMeans(T_null_c)
    T_null<-T_null_a+T_null_b+T_null_c
    ######################
    T_all_old<-cbind(T_0^2,T_null^2)
    kappa_old<-rowSums(T_null^2>T_0^2)
    tau_old<-apply(T_all_old,1,tau_cal)
    ord_old<-order(tau_old,decreasing = T)
    kappa_old_ord<-kappa_old[ord_old]
    
    if (v==1){
      cut_old<-min(which(kappa_old_ord!=0))
      selected_old<-rep(0,p)
      selected_old[ord_old[(1:p)<cut_old]]<-1
      
      FWER_v<-c(FWER_v,sum(selected_old*(beta_signal==0))>=1)
      
      PoWER_v<-c(PoWER_v,sum(selected_old*(beta_signal!=0))/sum((beta_signal!=0)))
      
      
      
    }else{
      
      selected_old<-rep(0,p)
      
      FWER_v<-c(FWER_v,sum(selected_old*(beta_signal==0))>=1)
      
      PoWER_v<-c(PoWER_v,sum(selected_old*(beta_signal!=0))/sum((beta_signal!=0)))
      
      
    }
    
    
    
    
    
 
    
    print(c(amplitude,iter))
    
 #   result.table.iner = rbind(result.table.iner,cbind(iter,FWER_old_a,PoWER_old_a,FWER_old,PoWER_old,FWER_v,PoWER_v))
 #   write.csv(result.table.iner,paste0('/oak/stanford/groups/zihuai/yuxinrui/zihuai/simulation_biometrics_summary/Summary_Feb_Power/simulation_result/amplitude_cs',amplitude,'_',p,'_rho_',rho,'.csv'),row.names=F,col.names=T,quote=F,sep='\t')
    
  }
  
  
  FWER_old<-c(FWER_old,mean(FWER_old_a))
  PoWER_old<-c(PoWER_old,mean(PoWER_old_a))

  FWER_old_sd<-c(FWER_old_sd,sd(FWER_old_a))
  PoWER_old_sd<-c(PoWER_old_sd,sd(PoWER_old_a))
  
  
  FWER_old_deran<-c(FWER_old_deran,mean(FWER_old_d))
  PoWER_old_deran<-c(PoWER_old_deran,mean(PoWER_old_d))
  
  FWER_old_deran_sd<-c(FWER_old_deran_sd,sd(FWER_old_d))
  PoWER_old_deran_sd<-c(PoWER_old_deran_sd,sd(PoWER_old_d))
  
  FWER_vanilla<-c(FWER_vanilla,mean(FWER_v))
  PoWER_vanilla<-c(PoWER_vanilla,mean(PoWER_v))
  
  FWER_vanilla_sd<-c(FWER_vanilla_sd,sd(FWER_v))
  PoWER_vanilla_sd<-c(PoWER_vanilla_sd,sd(PoWER_v))
  
  
  
    FWER_deran<-c(FWER_deran,mean(FWER_new_d))
  PoWER_deran<-c(PoWER_deran,mean(PoWER_new_d))
  
  FWER_deran_sd<-c(FWER_deran_sd,sd(FWER_new_d))
  PoWER_deran_sd<-c(PoWER_deran_sd,sd(PoWER_new_d))
  
  
  
  print(iter)
  
  
  
}


result.table = cbind(amp.seq,FWER_old_deran,PoWER_old_deran,FWER_deran,PoWER_deran,FWER_old,PoWER_old,FWER_vanilla,PoWER_vanilla,FWER_old_deran_sd,PoWER_old_deran_sd,FWER_deran_sd,PoWER_deran_sd,FWER_old_sd,PoWER_old_sd,FWER_vanilla_sd,PoWER_vanilla_sd)
data1 = result.table 



write.csv(result.table,paste0('/oak/stanford/groups/zihuai/yuxinrui/zihuai/simulation_biometrics_summary/Summary_Feb_Power/simulation_result/amplitude_new_cs_cderan_',p,'_rho_',rho,'amplitude_9.csv'),row.names=F,col.names=T,quote=F,sep='\t')






































