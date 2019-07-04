library(ggplot2)
library(reshape2)
library(GridLMM)
library(data.table)
library(readxl)
library(ggplot2)
library(dplyr)
library(lme4)
library(data.table)
library(stringr)

#I numericalized the genotypes, for the target region and the whole genotype table in tassle and ordered it the same as the phenotypes:
load('./Genotypes/GWASFilesHsftf9.Rimage')
#first update the v2 coordinates
z<-data.frame(fread('./Genotypes/V4FromV2.dblBed',fill=TRUE))
z$CHR<-sapply(z[,1],function(x) strsplit(x,' ')[[1]][1])
z$start<-sapply(z[,1],function(x) strsplit(x,' ')[[1]][3])
z$new<-z[,5]
z$SNP<-paste('S',z[,'CHR'],'_',z[,'start'],sep="")
gnK<-gnK[,which(colnames(gnK) %in% z$SNP)]

#Center Markers and make Kmat:
X<-gnK
X<-X*2
X<-as.matrix(X)
X[which(is.na(X))]<-0.5
#center it around 0
X<-X-1
X_centered = sweep(X,2,colMeans(X),'-') # center marker genotypes
#load in the high quality unimputed genos for K matrix
load('/Users/dangates/Desktop/HighFilteredReduced.Rimage')
genoMDS<-as.matrix(genoMDS)
genoMDS[which(is.na(genoMDS))]<-0.5
K_centered=sweep(genoMDS,2,colMeans(genoMDS),'-')
K = tcrossprod(K_centered) / ncol(K_centered)
rownames(K) = colnames(K) = rownames(genoMDS)
K = K/mean(diag(K))
#now I'm going to melt the DF to have different environments
meltb<-melt(data=data.frame(cbind(tb[,c('Taxa','YieldDrought','YieldHeat','YieldControl','YieldDroughtHeat')])))
meltc<-meltb[which(meltb$variable %in% c('YieldControl','YieldDrought')),]
colnames(meltc)<-c('Geno','variable','value')

#for ASI time:
#meltb<-melt(data=data.frame(cbind(tb[,c('Taxa','ASIDrought','ASIHeat','ASIControl','ASIDroughtHeat')])))
#meltc<-meltb[which(meltb$variable %in% c('ASIControl','ASIDrought')),]
#colnames(meltc)<-c('Geno','variable','value')

#for Flowering time:
#meltb<-melt(data=data.frame(cbind(tb[,c('Taxa','AnthesisDrought','AnthesisHeat','AntheisControl','AnthesisDroughtHeat')])))
#meltc<-meltb[which(meltb$variable %in% c('AntheisControl','AnthesisDrought')),]
#colnames(meltc)<-c('Geno','variable','value')


null_model = GridLMM_ML(formula = value~ 1 + variable +(1 + variable||Geno)+(1|Geno),data=meltc,initial_step = 0.01,mc.cores = 1)
h2_start = null_model$results[1,grep('.REML',colnames(null_model$results),fixed=T)]
V_setup = null_model$setup

greps<-paste('S',1:10,'_',sep="")
chunks<-lapply(greps,function(x){
  print(paste('Running',x))
  XXX<-X[,grep(x,colnames(X))]
  
  #now for the gxe_gwas
  gxe_gwas = GridLMM_GWAS(
    formula = value~ 1 + variable +(1 + variable||Geno)+(1|Geno),
    #formula = phenoMeanCentered~ 1 + envMeanCentered + year + tester + pop_1 + pop_2 + pop_3 + pop_4 + pop_5  + (1 | Geno/envMeanCentered) ,
    test_formula = ~1 + variable , # this is the model for each marker. ~1 means an intercept for the marker. ~1 + cov means an intercept plus a slope on `cov` for each marker
    reduced_formula = ~1, # This is the null model for each test. ~0 means just the error model. ~1 means the null includes the intercept of the marker, but not anything additional
    data = meltc, # The dataframe to look for terms from the 3 models
    weights = NULL, # optional observation-specific weights
    X = XXX, # The matrix of markers. Note: This can be of dimension n_g x p, where n_g is the number of unique genotypes.
    X_ID = 'Geno', # The column of the data that identifies each genotype. Each level of data$Geno should match a rowname of X
    h2_start = h2_start, # An optional vector of h2s to use as starting values for each random effect in the model. If NULL, will be calculated from the error model using GridLMM_ML
    h2_step = 1.1, # step size per random effect for testing alternate values of h2
    max_steps = 1, # maximum number of steps of size h2_step to take from h2_start
    relmat = list(Geno=K), # A list of Kernel matrices for the random effects. If X_ID (here Geno) is not included in this list, then it is calculated as tcrossprod(Xc)/ncol(Xc) where Xc is the centered (and optionally scaled) X. If any random effects are described in `error_model` but not provided here, the Kernel is assumed to be the identity matrix
    centerX = TRUE, # Should the markers be centered when calculating the GRM (only will be done if needed for calculating the GRM),
    scaleX = TRUE, # Should the markers be scaled to have constant variance when calculating the GRM?
    fillNAX = TRUE, # Should missing marker data be filled in with the mean allele frequency?
    V_setup = V_setup,
    method = 'REML', # Should the best model be selected by REML (if False, will be selected by ML)
    mc.cores = 1, # How many cores should be used for parallel processing. Unless X is large, tends to actually be faster with mc.cores = 1
    diagonalize = F,
    verbose = TRUE # Should progress be printed to the screen?
  )
})

SNPs<-unlist(sapply(chunks,function(x) return(x$results$p_value_REML.2)))
effectVs<-unlist(sapply(chunks,function(x) return(x$results$beta.4)))
#I need to double check that SNP names transfer like this:
names(SNPs)<-colnames(X)
names(effectVs)<-colnames(X)

colm<-data.frame(SNP=names(SNPs),CHR=sapply(names(SNPs),function(x) as.numeric(substr(strsplit(x,'_')[[1]][1],2,4))),BP=sapply(names(SNPs),function(x)as.numeric(strsplit(x,'_')[[1]][2])),P=SNPs,effect=effectVs)
colm<-colm[-which(is.na(colm$P)),]


#I also want to merge my Z and colm for downstream investigations
library(qqman)
merg<-merge(z[,c('SNP','new')],colm,by='SNP')
save(merg,file='DTMA_Anthesisgwas.Rimage')


