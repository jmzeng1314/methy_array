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

myLoad  
# 耗时步骤，运行一次后，就注释掉
if(F){
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
  dim(myNorm) 
  pD=myLoad$pd
  save(myNorm,pD,file = 'step2-champ_myNorm.Rdata')
}
load(file = 'step2-champ_myNorm.Rdata')
# 原来的450K经过质控过滤后是400K啦
beta.m=myNorm
group_list=myLoad$pd$Group
dim(beta.m) 
# 下面是表达矩阵标准3张图质量控制手段，生信技能树原创
if(T){
  
  
  dat=t(beta.m)
  dat[1:4,1:4] 
  library("FactoMineR")#画主成分分析图需要加载这两个包
  library("factoextra")  
  # 因为甲基化芯片是450K或者850K，几十万行的甲基化位点，所以PCA不会太快
  dat.pca <- PCA(dat , graph = FALSE) 
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               # palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  ggsave('all_samples_PCA.png')
  
  dat=beta.m
  dat[1:4,1:4] 
  cg=names(tail(sort(apply(dat,1,sd)),1000))#apply按行（'1'是按行取，'2'是按列取）取每一行的方差，从小到大排序，取最大的1000个
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #对那些提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵
  n=t(scale(t(dat[cg,]))) # 'scale'可以对log-ratio数值进行归一化
  n[n>2]=2 
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(group=group_list)
  rownames(ac)=colnames(n)  
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac,filename = 'heatmap_top1000_sd.png')
  dev.off()
  
  exprSet=beta.m
  pheatmap::pheatmap(cor(exprSet)) 
  # 组内的样本的相似性应该是要高于组间的！
  colD=data.frame(group_list=group_list)
  rownames(colD)=colnames(exprSet)
  pheatmap::pheatmap(cor(exprSet),
                     annotation_col = colD,
                     show_rownames = F,
                     filename = 'cor_all.png')
  dev.off() 
  exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
  dim(exprSet)
  # M=cor(log2(exprSet+1)) 
  M=cor(exprSet)
  pheatmap::pheatmap(M,annotation_col = colD)
  pheatmap::pheatmap(M,
                     show_rownames = F,
                     annotation_col = colD,
                     filename = 'cor_top500.png')
  dev.off() 
  
}


myLoad 
beta.m=myLoad$beta
# 下面是使用 wateRmelon 进行 归一化代码，目前被主流抛弃
if(F){
  # 使用 wateRmelon 进行 归一化
  library("wateRmelon")
  beta.m=beta.m[rowMeans(beta.m)>0.005,]
  pdf(file="rawBox.pdf")
  boxplot(beta.m,col = "blue",xaxt = "n",outline = F)
  dev.off()
  beta.m = betaqn(beta.m)
  pdf(file="normalBox.pdf")
  boxplot(beta.m,col = "red",xaxt = "n",outline = F)
  dev.off()
  
  # 然后进行简单的QC
  group_list=myLoad$pd$Group
  pdf(file="densityBeanPlot.pdf")
  par(oma=c(2,10,2,2))
  densityBeanPlot(beta.m, sampGroups = group_list)
  dev.off()
  pdf(file="mdsPlot.pdf")
  mdsPlot(beta.m, numPositions = 1000, sampGroups = group_list)
  dev.off()
  
  # 后续针对 beta.m 进行差异分析, 比如 minfi 包
  grset=makeGenomicRatioSetFromMatrix(beta.m,what="Beta")
  M = getM(grset)
  # 因为甲基化芯片是450K或者850K，几十万行的甲基化位点，统计检验通常很慢。
  dmp <- dmpFinder(M, pheno=group_list, type="categorical")
  dmpDiff=dmp[(dmp$qval<0.05) & (is.na(dmp$qval)==F),]
  dim(dmpDiff)
}



