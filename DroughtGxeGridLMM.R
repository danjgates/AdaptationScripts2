library(GridLMM)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

drought<-read.table('./GWASResiduals/Drought BLUPS-Table 1.csv',sep=",",header=TRUE)
phenotypes<-c(DA='dfm',PlantHeight='altpl',ASI='asi',DryWeight='pseco')
#bare cob weight, grain weight per hectare, weren't measured in 

#lapply(phenotypes,function(phenot){
phenot<-'pseco'

load('./Genotypes/HighFilteredReduced.Rimage')
genoMDS<-as.matrix(genoMDS)
genoMDS[which(is.na(genoMDS))]<-0.5
genoMDS<-genoMDS[which(rownames(genoMDS) %in% drought$Sample.ID),]
K_centered=sweep(genoMDS,2,colMeans(genoMDS,na.rm = TRUE),'-')
K = tcrossprod(K_centered) / ncol(K_centered)
rownames(K) = colnames(K) = rownames(genoMDS)
K = K/mean(diag(K))
dim(K)

drought<-drought[which(drought$Sample.ID %in% rownames(K)),]
#ok now pull out a phenotype from the table:
data<-subset(drought,Variable==phenot)
colnames(data)[which(colnames(data)=='Sample.ID')]<-'Geno'

null_model = GridLMM_ML(formula = BLUP~ 1 + Sitio + Manejo + (1 + Manejo || Geno) + (1|Geno),data=data,relmat = list(Geno=K),tolerance = 0.001,ML = F,REML=T,mc.cores=1)
h2_start = null_model$results[1,grep('.REML',colnames(null_model$results),fixed=T)]
V_setup = null_model$setup
h2_start

colGxEs<-lapply(1:10,function(chr){
  mat<-fread(paste('./Genotypes/Ch',chr,'Merged.hmp.txt',sep=""),skip=1,sep='\t',header=TRUE)
  mat<-data.frame(mat)
  row.names(mat)<-mat[,1]
  mat<-mat[,-1]
  mat<-mat[grep('SEED',rownames(mat)),]
  nms<-sapply(rownames(mat),function(x) strsplit(x,'.MR')[[1]][1])
  dups<-duplicated(nms)
  xx<-mat[-which(dups==TRUE),]
  rownames(xx)<-sapply(rownames(xx),function(x) strsplit(x,'.MR')[[1]][1])
  mat<-xx
  int<-intersect(rownames(mat),as.character(data$Geno))
  mat<-mat[int,]
  appls<-apply(mat,MARGIN=2,function(x) table(x)['0'])
  mat<-mat[,which(appls>25)] #rough filtering to reduce spurious stuff
  X<-mat
  X<-X[,-c(which(is.na(X),arr.ind=TRUE)[,2])] #I'm dropping any column w/ NA, since this is imputed it ends up being ~ 0.1% of columns
  X<-X*2
  X<-as.matrix(X)
  #X_centered = sweep(X,2,colMeans(X),'-') # center marker genotypes
  #cut into two: the first half and the second half
  brk<-as.integer(ncol(X)/5)
  chunks<-list(X[,1:brk],X[,(brk+1):(brk*2)],X[,(brk*2+1):(brk*3)],X[,(brk*3+1):(brk*4)],X[,(brk*4+1):ncol(X)])
  colGx<-lapply(chunks,function(X){
    #then iteratively hand it to this:
    map<-data.frame(snp=colnames(X),Chr=rep(chr,ncol(X)),pos=as.numeric(sapply(colnames(X),function(x) strsplit(x,'_')[[1]][2])))
    data<-data[data$Geno %in% rownames(X),]
    
    #note that the V_setup can be ~4GB so if you give that to a bunch of cores you're going to destroy your memory
    gxe_gwas = GridLMM_GWAS(
      formula = BLUP~ 1 + Sitio + Manejo + (1 + Manejo || Geno) + (1|Geno), # the same error model is used for each marker. It is specified similarly to lmer
      test_formula = ~1 + Sitio , # this is the model for each marker. ~1 means an intercept for the marker. ~1 + cov means an intercept plus a slope on `cov` for each marker
      reduced_formula = ~1, # This is the null model for each test. ~0 means just the error model. ~1 means the null includes the intercept of the marker, but not anything additional
      data = data, # The dataframe to look for terms from the 3 models
      weights = NULL, # optional observation-specific weights
      X = X, # The matrix of markers. Note: This can be of dimension n_g x p, where n_g is the number of unique genotypes.
      X_ID = 'Geno', # The column of the data that identifies each genotype. Each level of data$Geno should match a rowname of X
      h2_start = h2_start, # An optional vector of h2s to use as starting values for each random effect in the model. If NULL, will be calculated from the error model using GridLMM_ML
      h2_step = 0.1, # step size per random effect for testing alternate values of h2
      max_steps = 5, # maximum number of steps of size h2_step to take from h2_start
      X_map = map, # Optional. The marker positions.
      relmat = list(Geno=K), # A list of Kernel matrices for the random effects. If X_ID (here Geno) is not included in this list, then it is calculated as tcrossprod(Xc)/ncol(Xc) where Xc is the centered (and optionally scaled) X. If any random effects are described in `error_model` but not provided here, the Kernel is assumed to be the identity matrix
      centerX = TRUE, # Should the markers be centered when calculating the GRM (only will be done if needed for calculating the GRM),
      scaleX = TRUE, # Should the markers be scaled to have constant variance when calculating the GRM?
      fillNAX = FALSE, # Should missing marker data be filled in with the mean allele frequency?
      V_setup = V_setup,
      method = 'REML', # Should the best model be selected by REML (if False, will be selected by ML)
      mc.cores = 1, # How many cores should be used for parallel processing. Unless X is large, tends to actually be faster with mc.cores = 1
      diagonalize = F,
      verbose = TRUE # Should progress be printed to the screen?
    )
    
    pv<-gxe_gwas$results$p_value_REML.2
    colGxE<-data.frame(CHR=gxe_gwas$results$Chr,BP=gxe_gwas$results$pos,P=gxe_gwas$results$p_value_REML.2,SNP=gxe_gwas$results$snp,fitness=gxe_gwas$results$beta.5)
    return(colGxE)
  })
  
  library(qqman)
  colm<-rbind(colGx[[1]],colGx[[2]],colGx[[3]],colGx[[4]],colGx[[5]])
  return(colm)
})


colm<-colGxEs[[1]]
for( chr in 2:10 ){
  colm<-rbind(colm,colGxEs[[chr]])
}
colmD<-colm
trait<-names(phenotypes)[which(phenotypes==phenot)]
save(colmD,file=paste('DroughtGxE',trait,'.Rimage',sep=""))
#})
