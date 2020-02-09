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

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

load(file = 'anno_DEG.Rdata') 

gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] ) 
gene_down
gene_up
# 下面的函数依赖于kegg数据库连接网速，有时候中国大陆的朋友会失败
enrichKK <- enrichKEGG(gene         =  gene_up,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.1,
                      qvalueCutoff =0.1)
head(enrichKK)[,1:6] 
browseKEGG(enrichKK, 'hsa04512')
dotplot(enrichKK)
enrichKK=DOSE::setReadable(enrichKK, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
enrichKK 

#(3)可视化
#条带图
par(mfrow=c(2,1))
barplot(enrichKK,showCategory=20)
#气泡图
dotplot(enrichKK)
#下面的图需要映射颜色，设置和示例数据一样的geneList
geneList = deg$logFC
names(geneList)=deg$ENTREZID
geneList = sort(geneList,decreasing = T)
#(3)展示top5通路的共同基因，要放大看。
#Gene-Concept Network
cnetplot(enrichKK, categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
cnetplot(enrichKK, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
#Enrichment Map
emapplot(enrichKK)
#(4)展示通路关系,仅仅是针对于GO数据库结果。
# goplot(enrichKK)
#(5)Heatmap-like functional classification
heatplot(enrichKK,foldChange = geneList)


