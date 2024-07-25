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
Zscores_Ghost<-as.data.frame(fread("results_AD_Ghost_new2.csv"))
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
  chr_now<-Z_score_info$hg19.chr[index_k[1]]
  
  if(chr_now!=chr_old)
  {
    file_name<-paste0('Cor_',chr_now,".csv")
    Sigma_all<-as.matrix(fread(file_name))
    you_index<-Sigma_all[,1]
    Sigma<-Sigma_all[,-1]
    
    index_now<-which(Z_score_info$hg19.chr==chr_now)
    
    chr_old<-chr_now
  }
  
  aa<-which(index_now%in%index_k)
  Sigma_k<-as.matrix(Sigma[aa,aa],length(aa),length(aa))
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

######################
kappa_old<-apply(Z_all,1,kappa_1)
tau_old<-apply(Z_all,1,tau_1)

# tau_1(T_all_old[1,])

ord_old<-order(tau_old,decreasing = T)
kappa_old_ord<-kappa_old[ord_old]
print((cut_old<-min(which(kappa_old_ord!=0))))

selected_old<-rep(0,length(tau_old))
xuan<-ord_old[(1:length(tau_old))<cut_old]
selected_old[ord_old[(1:length(tau_old))<cut_old]]<-1

Z_select_info<-NULL
for(i in xuan)
{
  Z_select_info<-rbind(Z_select_info,Z_score_info[clusters==i,c(18,1:7,17)])
}
colnames(Z_select_info)[1]<-"cluster"
####################
Nameing<-read.csv("Select_name.csv")[,-1]
####################

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

xuan_less<-xuan[2]
for(i in 3:length(xuan))
{
  chr_ind<-(Z_pos[xuan_less,1]==Z_pos[xuan[i],1])
  pos_ind<-(abs(Z_pos[xuan_less,2]-Z_pos[xuan[i],2])<2e6)
  if(sum(chr_ind*pos_ind)==0)
  {
    xuan_less<-c(xuan_less,xuan[i])
  }
}
{
png(paste(c("Manhatton_Ghost_grp.png"),collapse = ""),1200,480)
par(mar=c(6,6,1,1))
plot(0,0,type="n",xlim=c(0,22),ylim=c(0,24),xlab = "",ylab="",axes=F)
mtext("chr",1,line=3,cex=2)
mtext(expression(tau[j]^'1/2'),2,line=3,cex=2)
axis(1,at=c(base_pos,22),labels = rep("",23))
axis(2,at=0:3*8,labels = 0:3*8)
abline(h=tau_old[ord_old[cut_old]]^(1/2),col="red",lty=2,lwd=3)
mtext(1:22,side=1,line=1,at=(base_pos+c(base_pos[-1],22))/2)
xy<-NULL
for(i in 1:length(tau_old))
{
  chr_index<-Z_pos[i,1]
  xx<-base_pos[chr_index]+LL[chr_index]*(Z_pos[i,2]-range_chr[chr_index,1])/(range_chr[chr_index,2]-range_chr[chr_index,1])
  points(xx,min(tau_old[i]^(1/2)),pch=20,col=c(gray(0.3),gray(0.6),"red","red")[selected_old[i]*2+1+chr_index%%2])
  xy<-rbind(xy,c(xx,min(tau_old[i]^(1/2))))
  if(i%%100==0)
  {print(i/100)}
}


library(shape)
i<-1
shuiping<-0
chuzhi<--2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,0+a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)


i<-2
shuiping<-0
chuzhi<--2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,0+a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-3
shuiping<-0
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-4
shuiping<-0
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-5
shuiping<-0
chuzhi<-7
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-6
shuiping<-0
chuzhi<-5
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-7
shuiping<-0
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-8
shuiping<--0.1
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-9
shuiping<-0
chuzhi<-4
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)


i<-10
shuiping<-0
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-11
shuiping<-0
chuzhi<-6
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-12
shuiping<-0.2
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)


i<-13
shuiping<-0
chuzhi<--1
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,0+a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-14
shuiping<-0
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-15
shuiping<--0.5
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-16
shuiping<-0
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-17
shuiping<-0.3
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)


i<-18
shuiping<-0.2
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-19
shuiping<-1.5
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-20
shuiping<-0
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-21
shuiping<-1.5
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-22
shuiping<-0
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-23
shuiping<-0.1
chuzhi<-1
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-24
shuiping<--0.3
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-25
shuiping<-0
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-26
shuiping<-1.5
chuzhi<-0.5
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-27
shuiping<-0
chuzhi<-2
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

i<-31
shuiping<-0
chuzhi<-5
You<-unique(Nameing$j[Nameing$cluster==xuan_less[i]])
for(a in 1:length(You))
{
  text(xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,You[a],adj=c(0.5,1-a*1.2))
}
Arrows(xy[xuan_less[i],1],xy[xuan_less[i],2],xy[xuan_less[i],1]+shuiping,xy[xuan_less[i],2]+chuzhi,arr.adj=1.2)

dev.off()
}

Z_select_info_final<-NULL
for(i in xuan)
{
  enen<-which(clusters==i)
  yy<-rep(NA,length(enen))
  yy[1]<-i
  hh<-rep(NA,length(enen))
  
  chi<-rep(NA,length(enen))
  chi[1]<-tau_old[i]
  
  You<-unique(Nameing$j[Nameing$cluster==i])
  You<-paste(You,collapse = ", ")
  hh[1]<-You
  hey<-paste(Z_score_info[enen,4],collapse = ", ")
  Z_select_info_final<-rbind(Z_select_info_final,c(i,Z_score_info[enen[1],3],hey,You,round(tau_old[i],2),NA))
}
write.table(Z_select_info_final,"Final.txt",sep="&",row.names = F)
write.csv(Z_select_info_final,"Final.csv",row.names = F)
####################
{
  Reference_info<-as.data.frame(fread("Reference_panel.csv"))
  
  reference_pos<-Reference_info[,c(4,6)]
  colnames(reference_pos)<-c("chr","pos")
  ss<-reference_pos$chr
  sss<-c()
  for(i in 1:length(ss))
  {
    ff<-strsplit(ss[i],split = "")[[1]]
    ddd<-paste(ff[-(1:3)],collapse = "")
    if(ddd%in%c("X","Y"))
    {
      ddd<-"23"
    }
    sss[i]<-as.integer(ddd)
  }
  reference_pos$chr<-sss
  reference_pos<-reference_pos[!is.na(sss),]
  Reference_info<-Reference_info[!is.na(sss),]
  
  
  
  Z_select_info_all<-NULL
  for(i in 1:nrow(Z_select_info))
  {
    aa<-which(Z_select_info$hg38.chr[i]==reference_pos$chr)
    # bb<-aa[which(Z_select_info$hg38.pos[i]==reference_pos$pos[aa])]
    kk<-Z_select_info$hg38.pos[i]-reference_pos$pos[aa]
    bb<-aa[which(abs(kk)==min(abs(kk)))]
    
    YY<-unique(Reference_info[bb,3])
    
    # Z_select_info_all_a<-NULL
    # Z_select_info_all_b<-NULL
    for(j in YY)
    {
      # Z_select_info_all_a<-rbind(Z_select_info_all_a,Z_select_info[i,])
      # Z_select_info_all_b<-rbind(Z_select_info_all_b,Reference_info[j,])
      Z_select_info_all<-rbind(Z_select_info_all,cbind(Z_select_info[i,],j))
    }
    
  }
  # colnames(Z_select_info)[1]<-"cluster"
  
}

####################
{
  Reference_info<-as.data.frame(fread("Reference_panel.csv"))
  
  reference_pos<-Reference_info[,c(14,15)]
  colnames(reference_pos)<-c("chr","pos")
  ss<-reference_pos$chr
  sss<-c()
  for(i in 1:length(ss))
  {
    ff<-strsplit(ss[i],split = "")[[1]]
    ddd<-paste(ff[-(1:3)],collapse = "")
    if(ddd%in%c("X","Y"))
    {
      ddd<-"23"
    }
    sss[i]<-as.integer(ddd)
  }
  reference_pos$chr<-sss
  reference_pos<-reference_pos[!is.na(sss),]
  Reference_info<-Reference_info[!is.na(sss),]
  
  
  
  Z_select_info_all<-NULL
  for(i in 1:nrow(Z_select_info))
  {
    aa<-which(Z_select_info$hg38.chr[i]==reference_pos$chr)
    # bb<-aa[which(Z_select_info$hg38.pos[i]==reference_pos$pos[aa])]
    kk<-Z_select_info$hg38.pos[i]-reference_pos$pos[aa]
    bb<-aa[which(abs(kk)==min(abs(kk)))]
    
    YY<-unique(Reference_info[bb,25])
    
    # Z_select_info_all_a<-NULL
    # Z_select_info_all_b<-NULL
    for(j in YY)
    {
      # Z_select_info_all_a<-rbind(Z_select_info_all_a,Z_select_info[i,])
      # Z_select_info_all_b<-rbind(Z_select_info_all_b,Reference_info[j,])
      Z_select_info_all<-rbind(Z_select_info_all,cbind(Z_select_info[i,],j))
    }
    
  }
  # colnames(Z_select_info)[1]<-"cluster"
  
}