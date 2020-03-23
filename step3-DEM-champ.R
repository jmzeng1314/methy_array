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
myLoad    # 存储了甲基化信号矩阵和表型信息。
load(file = 'step2-champ_myNorm.Rdata')
group_list=myLoad$pd$Group
table(group_list)
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
head(myDMP[[1]])
save(myDMP,file = 'step3-output-myDMP.Rdata')
# 还可以调动交互式界面修改阈值，调整差异化探针
# DMP.GUI(DMP=myDMP[[1]],beta=myNorm,group_list)

# 下面的分析, 非常的消耗计算资源
# 如果你有时间，就折腾
if(F){
  myDMR <- champ.DMR(beta = myNorm,pheno=group_list,method="Bumphunter")
  DMR.GUI(DMR=myDMR)
  
  myBlock <- champ.Block(beta = myNorm,pheno=group_list,arraytype="450K")
  head(myBlock$Block)
  Block.GUI(Block=myBlock,beta = myNorm,pheno=group_list,
            runDMP=TRUE,compare.group=NULL,arraytype="450K")
  
  myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]],
                       DMR=myDMR, arraytype="450K",adjPval=0.05, method="fisher")
  
  head(myGSEA$DMP)
  head(myGSEA$DMR)
  myEpiMod <- champ.EpiMod(beta=myNorm,pheno=group_list)
  
}




