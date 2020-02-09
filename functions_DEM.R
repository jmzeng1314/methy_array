visual_champ_DEM <- function(myLoad,myDMP,group_list,pro='test'){
  beta.m=myLoad$beta
  champDiff=myDMP[[1]] 
  head(champDiff)  
  colnames(champDiff)
  ## for volcano 
  if(T){
    nrDEG=champDiff
    head(nrDEG)
    attach(nrDEG)
    plot(logFC,-log10(P.Value))
    library(ggpubr)
    df=nrDEG
    df$v= -log10(P.Value) #df新增加一列'v',值为-log10(P.Value)
    ggscatter(df, x = "logFC", y = "v",size=0.5)
    
    df$g=ifelse(df$P.Value>0.05,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
                ifelse( df$logFC >0.5,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                        ifelse( df$logFC < -0.5,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
    )
    table(df$g)
    df$name=rownames(df)
    head(df)
    ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
    ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
              label = "name", repel = T,
              #label.select = rownames(df)[df$g != 'stable'] ,
              label.select =  head( rownames(df)[df$g != 'stable']),
              palette = c("#00AFBB", "#E7B800", "#FC4E07") )
    ggsave(paste0(pro,'_volcano.png'))
    
    ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
    df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                    ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
    table(df$p_c )
    ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
              palette = c("green", "red", "black") )
    ggsave(paste0(pro,'_MA.png'))
    
    
  }
  
  ## for heatmap 
  if(T){  
    dat=beta.m
    dat[1:4,1:4]
    table(group_list)
    deg=champDiff
    x=deg$logFC #deg取logFC这列并将其重新赋值给x
    names(x)=rownames(deg) #deg取probe_id这列，并将其作为名字给x
    cg=c(names(head(sort(x),100)),#对x进行从小到大排列，取前100及后100，并取其对应的探针名，作为向量赋值给cg
         names(tail(sort(x),100)))
    library(pheatmap)
    pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #对dat按照cg取行，所得到的矩阵来画热图
    n=t(scale(t(dat[cg,])))#通过“scale”对log-ratio数值进行归一化，现在的dat是行名为探针，列名为样本名，由于scale这个函数应用在不同组数据间存在差异时，需要行名为样本，因此需要用t(dat[cg,])来转换，最后再转换回来
    
    n[n>2]=2
    n[n< -2]= -2
    n[1:4,1:4]
    pheatmap(n,show_colnames =F,show_rownames = F)
    ac=data.frame(group=group_list)
    rownames(ac)=colnames(n) #将ac的行名也就分组信息（是‘no TNBC’还是‘TNBC’）给到n的列名，即热图中位于上方的分组信息
    pheatmap(n,show_colnames =F,
             show_rownames = F,
             cluster_cols = T, 
             annotation_col=ac,
             filename = paste0(pro,'_heatmap_top200_DEG_scale.png'))  
    pheatmap(dat[cg,],show_colnames =F,
             show_rownames = F,
             cluster_cols = T, 
             annotation_col=ac,
             filename = paste0(pro,'_heatmap_top200_DEG_raw.png'))  
    
  }
  
 
}