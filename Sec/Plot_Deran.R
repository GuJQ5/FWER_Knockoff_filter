library(data.table)
tau_1<-function(x)
{
  l<-sort(x,decreasing = T)
  return(l[1]-median(l[-1]))
}
kappa_1 <- function(x){
  l =x
  ls = which.max(l)-1
  
  return(ls)
}

# Zscores_Deran<-as.data.frame(fread('results_AD_Deran_16.txt'))
Zscores_Ghost<-as.data.frame(fread("results_AD_Deran_new2.csv"))
S_info<-as.data.frame(fread("Candidate_variants_info_other.csv"))
Zscores_Ghost_KO<-Zscores_Ghost[,-(1:18)]
Zscores_Ghost_Original<-Zscores_Ghost[,17]
Z_score_info<-Zscores_Ghost[,1:18]
Z_pos<-NULL

clusters<-Z_score_info$`clusters + cluster_count`
K_all<-length(unique(clusters))

Z_all<-NULL
all.rep.index<-NULL
chr_old<-0
rep_gene<-list()
for(k in 1:K_all)
{
  index_k<-which(clusters==k)
  chr_now<-Z_score_info$hg38.chr[index_k[1]]
  pos_now<-Z_score_info$hg38.pos[index_k]
  
  if(chr_now!=chr_old)
  {
    file_name<-paste0('Cor_',chr_now,".csv")
    Sigma_all<-as.matrix(fread(file_name))
    you_index<-Sigma_all[,1]
    Sigma<-Sigma_all[,-1]
    
    index_now<-which(Z_score_info$hg19.chr==chr_now)
    
    chr_old<-chr_now
    
  }
  info_index_now<-which(S_info$chr==chr_now)
  aa<-which(index_now%in%index_k)
  Sigma_k<-as.matrix(Sigma[aa,aa],length(aa),length(aa))
  # dd<-unlist(strsplit(S_info$HH_ref[index_now[aa]],split = ","))
  
  bb<-match(pos_now,S_info$pos[info_index_now])
  
  rep_gene[[k]]<-unique(unlist(strsplit(S_info$HH_ref[info_index_now[bb]],split = ",")))
  # rep.index<-which.max(rowSums(Sigma_k))
  # all.rep.index<-c(all.rep.index,index_k[rep.index])
  
  Z_0<-t(Zscores_Ghost_Original[index_k])%*%solve(Sigma_k)%*%Zscores_Ghost_Original[index_k]
  Z_KO<-rep(0,ncol(Zscores_Ghost_KO))
  for(m in 1:ncol(Zscores_Ghost_KO))
  {
    Z_KO[m]<-t(Zscores_Ghost_KO[index_k,m])%*%solve(Sigma_k)%*%Zscores_Ghost_KO[index_k,m]
  }
  Z_all<-rbind(Z_all,c(Z_0,Z_KO))
  Z_pos<-rbind(Z_pos,c(chr_now,mean(Z_score_info$hg19.pos[index_k])))
  print(k)
}

pi <- rep(0,nrow(Z_all))
gailv<- rep(0,nrow(Z_all))
for(i in 1:ncol(Zscores_Ghost_KO))
{
  T_all_old<-Z_all[,c(1,1+i)]
  kappa_old<-apply(T_all_old,1,kappa_1)
  tau_old<-apply(T_all_old,1,tau_1)
  
  ord_old<-order(tau_old,decreasing = T)
  kappa_old_ord<-kappa_old[ord_old]
  cut_old<-min(which(kappa_old_ord!=0))
  
  pi[ord_old[(1:length(tau_old))<cut_old]]<-pi[ord_old[(1:length(tau_old))<cut_old]]+1
  gailv[which(kappa_old_ord==0)]<-gailv[which(kappa_old_ord==0)]+1
  print(i)
}
eta<-0.98
pi<-pi/ncol(Zscores_Ghost_KO)
hist(pi)
selected_old_deran <- (pi>=eta)
xuan<-which(pi>=eta)
Z_select_info<-NULL
for(i in xuan)
{
  Z_select_info<-rbind(Z_select_info,Z_score_info[clusters==i,c(18,1:7,17)])
}
colnames(Z_select_info)[1]<-"cluster"
####################
Nameing<-read.csv("Select_name.csv")[,-1]
{
  range_chr<-NULL
  for(i in 1:22)
  {
    range_chr<-rbind(range_chr,range(Z_pos[Z_pos[,1]==i,2]))
  }
  LL<-range_chr[,2]-range_chr[,1]
  LL<-22*(LL/sum(LL))
  base_pos<-c(0,cumsum(LL)[-22])
}

xuan_less<-xuan[1:3]
# xuan_less[1]<-424

{
  png(paste(c("Manhatton_Deran_grp.png"),collapse = ""),1200,480)
  par(mar=c(6,6,1,1))
  plot(0,0,type="n",xlim=c(0,22),ylim=c(0,1),xlab = "",ylab="",axes=F)
  mtext("chr",1,line=3,cex=2)
  mtext("Selection Frequency",2,line=3,cex=2)
  axis(1,at=c(base_pos,22),labels = rep("",23))
  axis(2,at=0:5/5,labels = 0:5/5)
  abline(h=eta,col="red",lty=2,lwd=3)
  mtext(1:22,side=1,line=1,at=(base_pos+c(base_pos[-1],22))/2)
  xy<-NULL
  for(i in 1:length(tau_old))
  {
    chr_index<-Z_pos[i,1]
    xx<-base_pos[chr_index]+LL[chr_index]*(Z_pos[i,2]-range_chr[chr_index,1])/(range_chr[chr_index,2]-range_chr[chr_index,1])
    points(xx,pi[i],pch=20,col=c(gray(0.3),gray(0.6),"red","red")[selected_old_deran[i]*2+1+chr_index%%2])
    xy<-rbind(xy,c(xx,pi[i]))
    if(i%%100==0)
    {print(i/100)}
  }
  
  
  library(shape)
  i<-1
  shuiping<-0
  chuzhi<--0.1
  You<-"SELP"
  for(a in 1:length(You))
  {
    text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,0+a*1.2))
  }
  Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)
  
  
  i<-2
  shuiping<-0
  chuzhi<--0.1
  You<-c("TCF19", "PSORS1C3", "HCG27", "HLA-C")
  for(a in 1:length(You))
  {
    text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,0+a*1.2))
  }
  Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)
  
  i<-3
  shuiping<-0
  chuzhi<--0.1
  You<-"NXPE1"
  for(a in 1:length(You))
  {
    text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,0+a*1.2))
  }
  Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)
  
  i<-4
  shuiping<-0
  chuzhi<--0.1
  You<-c("APOE", "APOC1", "APOC1P1")
  for(a in 1:length(You))
  {
    text(xy[xuan[i],1]+shuiping,xy[xuan[i],2]+chuzhi,You[a],adj=c(0.5,0+a*1.2))
  }
  Arrows(xy[xuan[i],1],xy[xuan[i],2],xy[xuan[i],1]+shuiping,xy[xuan[i],2]+chuzhi,arr.adj=1.2)
  
  dev.off()
}
