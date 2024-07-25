

library(xtable)






pseq = c(50,100,200,500)

for(p in pseq){
  p=50
  nam <- paste("timeeff", p, sep = "")
  assign(nam, read.csv(paste0("time_efficient",p,".csv")))
 
}
Mean = NULL
for(p in pseq){

  name<-paste("m",p,sep="")
  assign(name,colMeans(eval(parse(text = paste0("timeeff",p)))))
  Mean = rbind(Mean,eval(parse(text = paste0("m",p))))
}
Mean.table = cbind(Mean[,2]+Mean[,3],Mean[,4],Mean[,5]+Mean[,6],Mean[,7])
table = format(round(Mean.table, 5), nsmall = 5) 
xtable(t(table))
sum.all = apply(Mean.table,1,sum)
xtable(sum.all)







for(p in pseq){

  nam <- paste("timeoriginal", p, sep = "")
  assign(nam, read.csv(paste0("time_original",p,".csv")))
  
}
Mean = NULL
for(p in pseq){
  
  name<-paste("moriginal",p,sep="")
  assign(name,colMeans(eval(parse(text = paste0("timeoriginal",p)))))
  Mean = rbind(Mean,eval(parse(text = paste0("moriginal",p))))
}
Mean.table = cbind(Mean[,2]+Mean[,3],Mean[,4],Mean[,5]+Mean[,6],Mean[,7])
table = format(round(Mean.table, 5), nsmall = 5) 
xtable(t(table))
sum.all = apply(Mean.table,1,sum)
print(sum.all)





for(p in pseq){
  
  nam <- paste("timederan", p, sep = "")
  assign(nam, read.csv(paste0("time_deran",p,".csv")))
  
}
Mean = NULL
for(p in pseq){
  
  name<-paste("mderan",p,sep="")
  assign(name,colMeans(eval(parse(text = paste0("timederan",p)))))
  Mean = rbind(Mean,eval(parse(text = paste0("mderan",p))))
}
Mean.table = cbind(Mean[,2]+Mean[,3],Mean[,4],Mean[,7],Mean[,8])
table = format(round(Mean.table, 5), nsmall = 5) 
xtable(t(table))
sum.all = apply(Mean.table,1,sum)
print(sum.all)




