# 一、TCGA-GBM数据下载（神经胶质瘤）
library(TCGAbiolinks)
library(SummarizedExperiment)
#查询TCGA中GBM的患者的RNA-seq基因表达数据，总共174个患者
query <- GDCquery(project = 'TCGA-GBM', 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
#下载上面搜索到的数据
GDCdownload(query, method = "api", files.per.chunk = 100)
#将下载查询的数据载入R
expdat <- GDCprepare(query = query)
#提取表达数据转换成数据框并保存
expr = assay(expdat)
expr = as.data.frame(expr)
save(expr,file = 'expr.Rdata')


# 二、提取mRNA数据
library(rtracklayer)
library(dplyr)
gtf <- import('Homo_sapiens.GRCh38.97.chr.gtf') 
gtf_df <- as.data.frame(gtf)  
#挑选ID、基因名、和基因功能
gene_df <- select(gtf_df,
                  c(gene_id,gene_name,gene_biotype))  
#去除重复
index <- duplicated(gene_d?f$gene_id) 
gene_df = gene_df[!index,]
dim(gene_df)
#挑选编码蛋白的基因
mRNA_df = gene_df[gene_df$gene_biotype == 'protein_coding',]  
dim(mRNA_df)
save(mRNA_df,file = 'mRNA_df.Rdata')
#载入之前的表达数据，根据上面提取的编码蛋白的mRNA的ID筛选表达数据中的mRNA
load("expr.Rdata")
exprSet = expr[match(mRNA_df$gene_id,rownames(expr)),]
dim(exprSet)
#剔除包含NA的行并保存
exprSet = na.omit(exprSet)
dim(exprSet)
save(exprSet, file = 'exprSet.Rdata')


# 三、ID转换，ensemblID转换为symbolID
library(limma)
load("exprSet.Rdata")
#用match方法，将id转换为mRNA矩阵中的基因名
exprSet$names = rownames(exprSet)
exprSet$names = mRNA_df[match(exprSet$names,mRNA_df$gene_id),2]
dim(exprSet)
#看是否有重复的基因名
table(duplicated(exprSet$names)) # 有3个重复的基因名
#用avereps相同ID取均值（limma包函数），对重复基因名取平均表达量，然后将基因名作为行名
#对数据除最后一列基因名（-ncol(exprSet)）之外的所有列的重复基因名取均值
exprSet = avereps(exprSet[,-ncol(exprSet)],ID = exprSet$names) 
dim(exprSet)
save(exprSet, file = 'exprSet_names_by_symbol.Rdata')


# 四、数据整理
#4.1 去除低表达量的基因
load("exprSet_names_by_symbol.Rdata")
#挑选一行中表达量等于0少于40的行，apply中1表示行，2表示列
pick_row <- apply(exprSet, 1, function(x){
  sum(x == 0) < 40
})
exprSet1 <- exprSet[pick_row,]
#4.2 分组(癌症组织和癌旁组织)
library(stringr)
#挑选出肿瘤和正常组患者barcode。tumor:01-09; normal:10-19; control:20-29
tumor <- colnames(exprSet1)[as.integer(substr(colnames(exprSet1),14,15)) < 10]
normal <- colnames(exprSet1)[as.integer(substr(colnames(exprSet1),14,15)) >= 10]
#将数据框分组
tumor_sample <- exprSet1[,tumor]
normal_sample <- exprSet1[,normal]
#然后再扩展列合并，后五个是正常
exprSet_by_group <- cbind(tumor_sample,normal_sample)
#添加分组标签
group_list <- c(rep('tumor',ncol(tumor_sample)),rep('normal',ncol(normal_sample)))
save(exprSet_by_group, group_list, file = 'exprSet_by_group_list.Rdata')


# 五、PCA
library(FactoMineR)
library(factoextra)
load("exprSet_by_group_list.Rdata")
#转置数据框变为pca用的形式，并添加一列分组因子
data = as.data.frame(t(exprSet_by_group))
data <- cbind(data,group = as.factor(group_list))
#pca
pca <- PCA(data[,-ncol(data)], graph = FALSE)
eig.val <- get_eigenvalue(pca)# 一列:特征值，二列:特征值的方差贡献度，三列:累计方差贡献度
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 30))
#pca作图
fviz_pca_ind(pca,
             geom.ind = "point", 
             col.ind = data$group, 
             palette = "jco", 
             addEllipses = TRUE, 
             legend.title = "Groups",
             title = 'PCA')


# 六、差异表达分析
#用limma和edgeR做差异分析
library(limma)
library(edgeR)
#edgeR
#创建DGEList类型变量
DGElist <- DGEList( counts = exprSet_by_group, group = factor(group_list))
#手动过滤,经验设置为cpm=1位为cutoff，也可以用filterByExpr自动过滤
keep_gene <- rowSums(cpm(DGElist) > 1 ) >= 2 
table(keep_gene)
DGElist <- DGElist[ keep_gene,keep.lib.sizes = FALSE]
#计算标准化因子
DGElist <- calcNormFactors(DGElist) # 计算归一化因子以对齐计数矩阵的列
DGElist$samples
#创建分组矩阵两列分别是正常，肿瘤。
design <- model.matrix( ~0 + factor(group_list))
rownames(design) <- colnames(DGElist)
colnames(design) <- levels(factor(group_list))
#limma
#转换RNA-Seq数据为线性建模做好准备
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
#给出一系列阵列，为每个基因拟合线性模型
fit <- lmFit(v, design)
#构造对应于一组参数的指定对比的对比矩阵。
cont.matrix <- makeContrasts(contrasts = c('tumor-normal'), levels = design)
#给定适合微阵列数据的线性模型，计算给定对比组的估计系数和标准误差
fit2 <- contrasts.fit(fit, cont.matrix)
#差分表达的经验贝叶斯统计
fit2 <- eBayes(fit2)
nrDEG_limma_voom = topTable(fit2, coef = 'tumor-normal', n = Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
save(nrDEG_limma_voom,file = 'nrDEG.Rdata')


# 七、火山图
library(ggplot2)
library(ggrepel)
load('nrDEG.Rdata')
nrDEG <- nrDEG_limma_voom
#设置阈值
nrDEG$change <- ifelse(nrDEG$adj.P.Val < 0.01 & abs(nrDEG$logFC) > 2.5,
                       ifelse(nrDEG$logFC > 2.5,'UP','DOWN'),
                       'NOT')
table(nrDEG$change)

#挑选重点关注基因，给这些基因在图上贴标签
nrDEG$sign <- ifelse(nrDEG$adj.P.Val < 0.0001 & abs(nrDEG$logFC) > 7,
                     rownames(nrDEG),
                     NA)
table(nrDEG$sign)
save(nrDEG,file='nrDEG_by_group.Rdata')
#画图
ggplot(data= nrDEG, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(alpha=0.8, size = 1) +
  theme_bw(base_size = 15) +
  theme(plot.title=element_text(hjust=0.5),   #  标题居中
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + # 网格线设置为空白
  geom_hline(yintercept=2 ,linetype=4) +
  geom_vline(xintercept=c(-2.5,2.5) ,linetype=4 ) +
  scale_color_manual(name = "", 
                     values = c("red", "green", "gray"),
                     limits = c("UP", "DOWN", "NOT")) +
  geom_label_repel(aes(label=sign), # 防止标签过多重叠
                   fontface="bold",
                   color="grey50",
                   box.padding=unit(0.35, "lines"),  # 文本框周边填充
                   point.padding=unit(0.5, "lines"), # 点周边填充
                   segment.colour = "grey50", # 连接点与标签的线段的颜色
                   force = T) + 
  labs(title = 'GBM DEG volcano')


# 八、热图
library("pheatmap")
nrDEG_Z = nrDEG[ order( nrDEG$logFC ), ]
nrDEG_F = nrDEG[ order( -nrDEG$logFC ), ]
#挑选logFC两头各50个基因作图
choose_gene = c( rownames( nrDEG_Z )[1:50], rownames( nrDEG_F )[1:50] )
#挑选所选的100个基因
choose_matrix = exprSet_by_group[ choose_gene, ]
#用scale函数将数据标准化，标准化是按行进行所以先转置
choose_matrix = t( scale( t( choose_matrix ) ) )
#设置+-2为最大阈值
choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2
#患者分组信息
annotation_col = data.frame( group = factor( group_list ) )
rownames( annotation_col ) = colnames( exprSet_by_group )
#画图
pheatmap( fontsize_row = 4,
          choose_matrix, 
          annotation_col = annotation_col,  # 患者分组信息作为列注释
          show_rownames = T,
          show_colnames = F,
          annotation_legend = T, 
          cluster_cols = T,
          filename = 'heatmap.png')


# 七、富集分析
#1。kegg富集上调下调基因
#转换ID
library( "clusterProfiler" )
library( "org.Hs.eg.db" )
load("nrDEG_by_group.Rdata")
nrDEG$SYMBOL <- rownames(nrDEG)
df <- bitr( rownames( nrDEG ), fromType = "SYMBOL", toType = c( "ENTREZID" ), 
            OrgDb = org.Hs.eg.db )
head( df )
nrDEG = merge( nrDEG, df, by = 'SYMBOL' )
head( nrDEG )
#挑选上调下调的差异表达基因
gene_up = nrDEG[ nrDEG$change == 'UP', 'ENTREZID' ] 
gene_down = nrDEG[ nrDEG$change == 'DOWN', 'ENTREZID' ]
gene_diff = c( gene_up, gene_down )
gene_all = as.character(nrDEG[ ,'ENTREZID'] )
g_list = list( gene_up = gene_up, gene_down = gene_down, gene_diff = gene_diff)
#用clusterProfiler包的KEGG在线注释
kk.up <- enrichKEGG(gene = gene_up,
                    organism = 'hsa',
                    universe = gene_all,
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01)
kk.dowm <- enrichKEGG(gene = gene_down,
                      organism = 'hsa',
                      universe = gene_all,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01)
#将KEGG富集注释结果转化为数据框
kegg_down_dt <- as.data.frame(kk.dowm)
kegg_up_dt <- as.data.frame( kk.up )
#以p<0.05的阈值挑选显著通路，加上下调分组信息并合并
down_kegg <- kegg_down_dt[ kegg_down_dt$pvalue < 0.05, ]
down_kegg$group <- 'down_pathway'
up_kegg <- kegg_up_dt[ kegg_up_dt$pvalue < 0.05, ]
up_kegg$group <- 'up_pathway'
dat = rbind(up_kegg,down_kegg)
#对p值取负对数，将分组因子化
dat$pvalue = -log10(dat$pvalue)
dat$group =  factor(dat$group)
#用ggpubr画简单柱状图
library(ggpubr)
ggbarplot(dat,x = 'Description',y = 'pvalue',
          fill = 'group',
          color = 'white',
          palette = 'jco',
          sort.val = 'asc',
          xlab = 'Pathway names',
          ylab = '-log10 P-value',
          title = 'Pathway enrichment') +
  rotate() +
  theme_minimal()

#2。go富集差异基因
#用clusterProfiler包的go在线注释
#
BP <- enrichGO( gene          =  gene_diff,
                universe      =  gene_all,
                OrgDb         =  org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           =  'BP',
                pAdjustMethod = "BH",
                pvalueCutoff  =  0.01,
                qvalueCutoff  =  0.01,
                readable      =  TRUE)
barplot(BP,showCategory=20)
dotplot(BP,showCategory=20)
ALL <- enrichGO( gene          =  gene_diff,
                universe      =  gene_all,
                OrgDb         =  org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           =  'ALL',
                pAdjustMethod = "BH",
                pvalueCutoff  =  0.01,
                qvalueCutoff  =  0.01,
                readable      =  TRUE)
