
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

load(file = 'step3-output-myDMP.Rdata')
deg=myDMP[[1]]
head(deg)
length(unique(deg$gene)) 
deg$g=ifelse(abs(deg$logFC) < 0.2,'stable',
             ifelse(deg$logFC > 0.2,'UP','DOWN'))
table(deg$g)
head(deg)
deg$symbol=deg$gene
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)

DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
save(DEG,file = 'anno_DEG.Rdata')


gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] ) 
source('kegg_and_go_up_and_down.R')
# 有时候kegg数据库抽风，这个函数会失败，这个锅，Y叔背
run_kegg(gene_up,gene_down,pro='test_methy')
# 需要多go数据库的3个条目进行3次富集分析，非常耗时。
# 所以我注释掉了下面的代码
# run_go(gene_up,gene_down,pro='npc_VS_normal')

# 下面的GO肯定会慢，但是会成功运行
go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
library(ggplot2)
library(stringr)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave('gene_up_GO_all_barplot.png') 
go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave('gene_down_GO_all_barplot.png')


