## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2020-02-09 16:46:35
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2020-02-09   First version
###
### ---------------



rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-output.Rdata')
beta.m=myLoad$beta
# 后续针对 beta.m 进行差异分析, 比如 minfi 包
grset=makeGenomicRatioSetFromMatrix(beta.m,what="Beta")
M = getM(grset)
group_list=myLoad$pd$Group
# 因为甲基化芯片是450K或者850K，几十万行的甲基化位点，统计检验通常很慢。
dmp <- dmpFinder(M, pheno=group_list, type="categorical")
dmpDiff=dmp[(dmp$qval<0.05) & (is.na(dmp$qval)==F),]
dim(dmpDiff)

load(file = 'step3-output-myDMP.Rdata')
champDiff=myDMP[[1]]

dim(dmpDiff)
dim(champDiff)
length(intersect(rownames(dmpDiff),rownames(champDiff)))

source('functions_DEM.R')
visual_champ_DEM(myLoad,myDMP,group_list,pro='test')



