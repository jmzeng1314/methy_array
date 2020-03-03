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

require(GEOquery)
require(Biobase)
library("impute")

# 这里举例的  group.txt 和 data.txt 是自己截图的6个甲基化芯片数据，2个分组
# 方便走差异分析流程。

# 任何一个GEO的数据集，都可以自行这里  group.txt 和 data.txt 文件
# 或者走后面的 GEOquery 流程，取决于你自己的需求哈

info=read.table("group.txt",sep="\t",header=T)
library(data.table)
b=info
rownames(b)=b[,1]
# 如果你的甲基化信号矩阵，自己在Excel表格里面整理好。
# 就走下面的fread流程
a=fread("data.txt",data.table = F )
a[1:4,1:4]
rownames(a)=a[,1]
a=a[,-1]
beta=as.matrix(a)
beta=impute.knn(beta)
betaData=beta$data
betaData=betaData+0.00001
a=betaData
a[1:4,1:4]
identical(colnames(a),rownames(b))
# 一定要保证，甲基化信号值矩阵，和表型信息，是一一对应的

library(ChAMP)
# beta 信号值矩阵里面不能有NA值
myLoad=champ.filter(beta = a,pd = b)
myLoad
save(myLoad,file = 'step1-output.Rdata')

# 如果你使用GEO数据库下载甲基化信号值矩阵文件
# 下面的代码你也需要理解哦。
if(F){
  require(GEOquery)
  require(Biobase)
  eset <- getGEO("GSE68777",destdir = './',AnnotGPL = T,getGPL = F)
  beta.m <- exprs(eset[[1]])
  ## 顺便把临床信息制作一下，下面的代码，具体每一个项目都是需要修改的哦
  pD.all <- pData(eset[[1]])
  pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
  head(pD)
  names(pD)[c(3,4)] <- c("group", "sex")
  pD$group <- sub("^diagnosis: ", "", pD$group)
  pD$sex <- sub("^Sex: ", "", pD$sex)

  library(ChAMP)
  # beta 信号值矩阵里面不能有NA值
  myLoad=champ.filter(beta = beta.m ,pd = pD)
  myLoad
  save(myLoad,file = 'step1-output.Rdata')
}

# 两种方法，都是为了制作 champ 的对象
# 后续分析，使用这个myLoad变量即可
