#source("GhostKnockoff.R")

##### Premelianry, including source and library

#Enable command line arguments


#amplitude_value = as.numeric(args[2])

#args = commandArgs(TRUE)
#p = as.numeric(args[1])
#rho = as.numeric(args[2])

p =500

rho =0.5


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
amp.seq = c(8)





#Sigma<-toeplitz(rho^(1:p-1))
#diag(Sigma)<-1

###### Compound Symmetry
Sigma<-diag(1-rho,p)+rho


Sigma_chol<-chol(Sigma)
SigmaInv<-solve(Sigma)

s.seq = list()
PP.seq = list()
V1_seq = list()
V1.right.seq = list()
V2.right.seq = list()
for(i in 1:42){
  s.seq[[i]] = list()
  PP.seq[[i]] = list()
  V1_seq[[i]] = list()
  V1.right.seq[[i]] = list()
  V2.right.seq[[i]] = list()
  
}

for(i in 1:42){
  
  M.19<-18+i
  s<-create.solve_sdp_M(Sigma,M=M.19)
  D<-Matrix(diag(s,length(s)))
  PP.19<-(diag(1,p)-D%*%SigmaInv)
  V1.19<-(M.19+1)/M.19*D-D%*%SigmaInv%*%D
  V1.19<-(V1.19+t(V1.19))/2
  V1.right.19<-chol(V1.19)
  V2.right.19<-chol(D)
  
  s.seq[[i]] = s
  PP.seq[[i]] = PP.19
  V1_seq[[i]] = V1.19
  V1.right.seq[[i]] = V1.right.19
  V2.right.seq[[i]] = V2.right.19
  
  
  
}

v.seq = c(rep(1,20),rep(2,20),rep(3,2))


FWER= list()
Power = list()
for(i in 1:42){
  FWER[[i]] = list()
  Power[[i]] = list()
  
}

for (amplitude in amp.seq){
  
  FWER_M = list()
  Power_M = list()
  for (i in 1:42){
    FWER_M[[i]] = list()
    Power_M[[i]] = list()
    
  }
  
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
    
    
    
    for(i in 1:42){
      M.19 = 18+i
      M = 18+i
  
      PP.19 = PP.seq[[i]] 
      V1.19 = V1_seq[[i]] 
      V1.right.19 = V1.right.seq[[i]] 
      V2.right.19 = V2.right.seq[[i]] 
      v = v.seq[i]
      T_null_a<-(PP.19%*%T_0)%x%t(rep(1,M.19))
      T_null_b<-(t(V1.right.19)%*%rnorm(p))%x%t(rep(1,M))
      T_null_c<-t(V2.right.19)%*%matrix(rnorm(p*M),p,M)
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
      
      FWER_M[[i]] = c(FWER_M[[i]] ,sum(selected_old*(beta_signal==0))>=1)
      Power_M[[i]] =c(Power_M[[i]],sum(selected_old*(beta_signal!=0))/sum((beta_signal!=0)))
      
    }
    
    
    
    print(c(amplitude,iter))
    
    # result.table.iner = rbind(result.table.iner,cbind(iter,FWER_new.19,PoWER_new.19,FWER_new.24,PoWER_new.24,FWER_new.27,PoWER_new.27,FWER_new.29,PoWER_new.29,FWER_new.31,PoWER_new.31,FWER_new.34,PoWER_new.34,FWER_new.37,PoWER_new.37,FWER_new.39,PoWER_new.39,FWER_new.41,PoWER_new.41,FWER_new.44,PoWER_new.44,FWER_new.47,PoWER_new.47,FWER_new.49,PoWER_new.49))
    
    #  write.csv(result.table.iner,paste0('/oak/stanford/groups/zihuai/yuxinrui/zihuai/simulation_biometrics_summary/Summary_Feb_Power/simulation_M/amplitude_nv',amplitude,'_',p,'_rho_',rho,'.csv'),row.names=F,col.names=T,quote=F,sep='\t')
    
  }
  
  
  for(i in 1:42){
    FWER[[i]]<-c(FWER[[i]],mean(as.numeric(FWER_M[[i]])))
    Power[[i]]<-c(Power[[i]],mean(as.numeric(Power_M[[i]])))
    
  }
 
  

  print(iter)
  
  
  
}


#result.table = cbind(amp.seq,FWER.19,PoWER.19,FWER.22,PoWER.22,FWER.24,PoWER.24,FWER.26,PoWER.26,FWER.29,PoWER.29,FWER.31,PoWER.31,FWER.34,PoWER.34,FWER.37,PoWER.37,FWER.39,PoWER.39,FWER.41,PoWER.41,FWER.44,PoWER.44,FWER.47,PoWER.47,FWER.49,PoWER.49)
#data1 = result.table 



write.csv(FWER,paste0('/oak/stanford/groups/zihuai/yuxinrui/zihuai/simulation_biometrics_summary/Summary_Feb_Power/simulation_M/FWER_cs_nv_',p,'_rho_',rho,'.csv'),row.names=F,col.names=T,quote=F,sep='\t')


write.csv(Power,paste0('/oak/stanford/groups/zihuai/yuxinrui/zihuai/simulation_biometrics_summary/Summary_Feb_Power/simulation_M/Power_cs_nv_',p,'_rho_',rho,'.csv'),row.names=F,col.names=T,quote=F,sep='\t')






































