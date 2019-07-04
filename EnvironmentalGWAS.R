library(data.table)


#Pull in the Environmental matrix
load('./GWASResiduals/WorldClimMat.Rimage')

#specify the layers you want to run (an example from a rerun is hashed out below)
#layers<-c("cld6MonthGrowAvg","dtr6MonthGrowAvg","frs6MonthGrowAvg","pet6MonthGrowAvg","pH","pre6MonthGrowAvg","tmn6MonthGrowAvg","tmp6MonthGrowAvg","tmx6MonthGrowAvg","vap6MonthGrowAvg","Waterlog.Percent","wet6MonthGrowAvg")
layer<-'cld6MonthGrowAvg'
#make a loop to run through the manhattans:
#sapply(layers,function(layer){
#layer<-'altitude'
precip<-finalMat[,layer]
names(precip)<-finalMat$SampleID
#build up some environmental GWAS

#pop structure:
load('./Genotypes/HighFilteredReduced.Rimage')
inters<-intersect(names(precip),rownames(genoMDS))
genoMDS<-genoMDS[inters,]
precip<-precip[inters]
pops<-cmdscale(dist(genoMDS),k=5)


colms<-lapply(1:10,function(chr){
mat<-data.frame(fread(paste('./Genotypes/Ch',chr,'Merged.hmp.txt',sep="")))
nms<-sapply(mat[,1],function(x) strsplit(x,'.M')[[1]][1])
dups<-duplicated(nms)
xx<-mat[-which(dups==TRUE),]
rownames(xx)<-sapply(xx[,1],function(x) strsplit(x,'.M')[[1]][1])
mat<-xx
mat<-mat[,-1]
mat<-mat[names(precip),]
mono<-apply(mat,MARGIN=2,function(x) length(table(x)))
mat<-mat[,-which(mono==1)]

#ok read in SNPs:
colm<-sapply(1:ncol(mat),function(x){
mergetab<-data.frame(cbind(geno=mat[,x],Precip=precip,V1=pops[,1],V2=pops[,2],V3=pops[,3],V4=pops[,4],V5=pops[,5]))
pv<-summary(lm(Precip~geno+V1+V2+V2+V4+V5,data=mergetab))$coefficient[2,4]
pred<-predict(lm(Precip~geno+V1+V2+V2+V4+V5,data=mergetab),data.frame(geno=1,V1=0,V2=0,V3=0,V4=0,V5=0))-predict(lm(Precip~geno+V1+V2+V3+V4+V5,data=mergetab),data.frame(geno=0,V1=0,V2=0,V3=0,V4=0,V5=0))
return(c(CHR=chr,BP=as.numeric(strsplit(colnames(mat)[x],'_')[[1]][2]),P=pv,prediction=pred))
})
return(colm)
})
colm<-data.frame(rbind(t(colms[[1]]),t(colms[[2]]),t(colms[[3]]),t(colms[[4]]),t(colms[[5]]),t(colms[[6]]),t(colms[[7]]),t(colms[[8]]),t(colms[[9]]),t(colms[[10]])))
library(dplyr)
colm<-mutate(colm,SNP=paste('S',1:nrow(colm),sep="_"))
library(qqman)
save(colm,file=paste(layer,'GWAS_2.Rimage',sep=""))
#manhattan(colm)
#})



