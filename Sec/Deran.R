library(data.table)
library(Matrix)
library(susieR)
# library(gdsfmt)
source('GhostKnockoff_hello.R')

#load z-scores
M<-1 # number of knockoffs
M.deran<-50
Zscores<-as.data.frame(fread('AD_Zscores_candidate.csv'))
# Zscores<-as.data.frame(fread('AD_Zscores_candidate.csv'))
Zscores_all<-NULL
#options(warn=2)
cluster_count<-0
for(chr in 1:22)
{
  # chr<-1
  file_name<-paste0('Cor_',chr,".csv")
  Sigma_all<-as.matrix(fread(file_name))
  you_index<-Sigma_all[,1]
  Sigma<-Sigma_all[,-1]
  
  Sigma.distance = as.dist(1 - abs(Sigma))
  fit = hclust(Sigma.distance, method="single")
  corr_max = 0.75
  clusters = cutree(fit, h=1-corr_max)
  K_all<-length(unique(clusters))
  # cluster.index<-c()
  # for(k in 1:K_all)
  # {
  #   temp.cor<-Sigma[clusters==k,clusters==k,drop=F]
  #   if(sum(clusters==k)==1){
  #     cluster.index<-c(cluster.index,which(clusters==k))
  #   }else{
  #     diag(temp.cor)<-NA
  #     cluster.index<-c(cluster.index,which(clusters==k)[which.max(apply(abs(temp.cor),1,mean,na.rm=T))])
  #   }
  # }
  
  
  Sigma<-forceSymmetric(Sigma)
  SigmaInv<-solve(Sigma)
  Zscores_chr<-Zscores[Zscores$hg38.chr==chr,]
  Zscores_chr<-Zscores_chr[you_index,]
  # Zscores_chr<-Zscores_chr[cluster.index,]
  Z_chr<-Zscores_chr$Z
  
  s<-create.solve_group_sdp_pca_M(as.matrix(Sigma),M=M,clusters=clusters)
  diag_s<-s$S
  # diag_s<-Matrix(diag(s,length(s)))
  P.each<-diag(1,ncol(Sigma))-diag_s%*%SigmaInv
  Sigma_k<-2*diag_s - diag_s%*%SigmaInv%*%diag_s
  #new trick by Jiaqi to improve computing time
  L1.left<-t(chol(forceSymmetric(Sigma_k-(M-1)/M*diag_s)))
  Normal_1<-L1.left%*%matrix(rnorm(M.deran*ncol(L1.left)),ncol(L1.left),M.deran)
  
  GK.Zscore_0<-as.matrix(Z_chr,ncol=1)
  GK.Zscores_k<-as.vector(P.each%*%GK.Zscore_0)+Normal_1
  colnames(GK.Zscores_k)<-paste0('Z_',1:ncol(GK.Zscores_k))
  
  temp.results<-cbind(Zscores_chr,clusters+cluster_count,as.matrix(GK.Zscores_k))
  Zscores_all<-rbind(Zscores_all,temp.results)
  cluster_count<-cluster_count+K_all
  print(chr)
}
# out.dir<-'/oak/stanford/groups/zihuai/jiaqigu/AD2/'

# temp.filename<-paste0(out.dir,'results_AD_Ghost_new_',replicate,'.txt')
write.csv(Zscores_all,'results_AD_Deran_new2.csv',row.names=F)

print('done!')
