# 待修改，增加背景基因集
enrichment <- function(geneset, ont='MF'){
  up_overlap <- geneset
  up <- bitr( up_overlap, fromType = "SYMBOL", toType = c( "ENTREZID" ), 
              OrgDb = org.Hs.eg.db )
  kegg <- enrichKEGG(gene          = up$ENTREZID,
                      keyType       = "kegg",
                      organism      = 'hsa',
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05)
  go <- enrichGO( gene         = up$ENTREZID,
                     OrgDb         =  org.Hs.eg.db,
                     keyType       = 'ENTREZID',
                     ont           =  ont,
                     pAdjustMethod = "BH",
                     pvalueCutoff  =  0.05,
                     qvalueCutoff  =  0.05,
                     readable      =  TRUE)
  return(list(kegg=kegg,go=go))
}


number_ticks <- function(n) {function(limits) pretty(limits, n)}


monk.replicationPlot <- function(humandata, monkeydata, xlabel, ylabel, xlimit, ylimit){
  homologs <- intersect(humandata$genename, monkeydata$genename)
  length(homologs)
  human_homolog <- humandata[humandata$genename %in% homologs,]
  monkey_homolog <- monkeydata[monkeydata$genename %in% homologs,]
  # 去重保留第一次出现的基因名（本质上是随机保留相同基因
  monkey_homolog <- monkey_homolog[!duplicated(monkey_homolog$genename),]
  homo_data <- inner_join(human_homolog, monkey_homolog, by="genename")
  homo_data$change <- ifelse(homo_data$change.x=="UP" & homo_data$change.y=="UP", "BOTH", 
                             ifelse(homo_data$change.x=="UP", "HUMAN",
                                    ifelse(homo_data$change.y=="UP", "MONKEY", "NOT")))
  ahba.array <- homo_data$logFC 
  rep.array  <- homo_data$log2FoldChange
  print(cor.test(ahba.array, rep.array))
  ahba.sig.genes <- homo_data[homo_data$change=="HUMAN",]$genename
  rep.sig.genes <- homo_data[homo_data$change=="MONKEY",]$genename
  
  # Color significant AHBA genes in red
  color.arr        <- rep('grey', length(ahba.array))
  color.arr[which(homo_data$change=="BOTH")] <- 'springgreen4'
  color.arr[which(homo_data$change=="MONKEY")]   <- 'gold'
  color.arr[which(homo_data$change=="HUMAN")]   <- 'deepskyblue2'
  
  plot_df <- as.data.frame(cbind(as.numeric(ahba.array),  as.numeric(rep.array), color.arr))
  colnames(plot_df) <- c('AHBA','RepDat', 'color')
  # Keep the proper color plotting order
  plot_df$color <- factor(plot_df$color, levels = c('grey','gold','springgreen4','deepskyblue2'))
  plot_df$AHBA   <- as.numeric(as.character(plot_df$AHBA))
  plot_df$RepDat <- as.numeric(as.character(plot_df$RepDat))
  for (col in levels(plot_df$color) ){
    field.ahba   <- paste(col, 'AHBA', sep = '')
    field.repdat <- paste(col, 'RepDat', sep = '')
    
    plot_df[[field.ahba]]   <- plot_df$AHBA
    plot_df[[field.repdat]] <- plot_df$RepDat
    plot_df[[field.ahba]][which(!plot_df$color == col)] <- NA
    plot_df[[field.repdat]][which(!plot_df$color == col)] <- NA
  }
  ols <- lm(RepDat ~ AHBA,
            data = plot_df)
  ggplot( plot_df ) +
    geom_point( aes(x = greyAHBA, y = greyRepDat, color=color), size = 2.5) +
    geom_point( aes(x = goldAHBA, y = goldRepDat, color=color), size = 2.5) +
    geom_point( aes( x = deepskyblue2AHBA, y = deepskyblue2RepDat, color = color), size = 2.5) + 
    geom_point( aes(x = springgreen4AHBA, y = springgreen4RepDat, color = color), size = 2.5) +
    xlab(xlabel) + ylab(ylabel) +
    scale_color_manual(breaks=unique(plot_df$color), values=as.character(unique(plot_df$color))) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_x_continuous(limits = c(-1*xlimit, xlimit), breaks = number_ticks(11))  +
    scale_y_continuous(limits = c(-1*ylimit, ylimit), breaks = number_ticks(11))  +
    geom_abline(intercept = ols$coefficients[1],
                slope = ols$coefficients[2])
}


hum.replicationPlot <- function(humandata, monkeydata, xlabel, ylabel, xlimit, ylimit){
  homologs <- intersect(humandata$genename, monkeydata$genename)
  length(homologs)
  human_homolog <- humandata[humandata$genename %in% homologs,]
  monkey_homolog <- monkeydata[monkeydata$genename %in% homologs,]
  # 去重保留第一次出现的基因名（本质上是随机保留相同基因
  monkey_homolog <- monkey_homolog[!duplicated(monkey_homolog$genename),]
  homo_data <- inner_join(human_homolog, monkey_homolog, by="genename")
  homo_data$change <- ifelse(homo_data$change.x=="UP" & homo_data$change.y=="UP", "BOTH", 
                             ifelse(homo_data$change.x=="UP", "HUMAN",
                                    ifelse(homo_data$change.y=="UP", "MONKEY", "NOT")))
  ahba.array <- homo_data$logFC.x 
  rep.array  <- homo_data$logFC.y
  print(cor.test(ahba.array, rep.array))
  ahba.sig.genes <- homo_data[homo_data$change=="HUMAN",]$genename
  rep.sig.genes <- homo_data[homo_data$change=="MONKEY",]$genename
  
  # Color significant AHBA genes in red
  color.arr        <- rep('grey', length(ahba.array))
  color.arr[which(homo_data$change=="BOTH")] <- 'springgreen4'
  color.arr[which(homo_data$change=="MONKEY")]   <- 'gold'
  color.arr[which(homo_data$change=="HUMAN")]   <- 'deepskyblue2'
  
  plot_df <- as.data.frame(cbind(as.numeric(ahba.array),  as.numeric(rep.array), color.arr))
  colnames(plot_df) <- c('AHBA','RepDat', 'color')
  # Keep the proper color plotting order
  plot_df$color <- factor(plot_df$color, levels = c('grey','gold','springgreen4','deepskyblue2'))
  plot_df$AHBA   <- as.numeric(as.character(plot_df$AHBA))
  plot_df$RepDat <- as.numeric(as.character(plot_df$RepDat))
  for (col in levels(plot_df$color) ){
    field.ahba   <- paste(col, 'AHBA', sep = '')
    field.repdat <- paste(col, 'RepDat', sep = '')
    
    plot_df[[field.ahba]]   <- plot_df$AHBA
    plot_df[[field.repdat]] <- plot_df$RepDat
    plot_df[[field.ahba]][which(!plot_df$color == col)] <- NA
    plot_df[[field.repdat]][which(!plot_df$color == col)] <- NA
  }
  ols <- lm(RepDat ~ AHBA,
            data = plot_df)
  ggplot( plot_df ) +
    geom_point( aes(x = greyAHBA, y = greyRepDat, color=color), size = 2.5) +
    geom_point( aes(x = goldAHBA, y = goldRepDat, color=color), size = 2.5) +
    geom_point( aes( x = deepskyblue2AHBA, y = deepskyblue2RepDat, color = color), size = 2.5) + 
    geom_point( aes(x = springgreen4AHBA, y = springgreen4RepDat, color = color), size = 2.5) +
    xlab(xlabel) + ylab(ylabel) +
    scale_color_manual(breaks=unique(plot_df$color), values=as.character(unique(plot_df$color))) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    scale_x_continuous(limits = c(-1*xlimit, xlimit), breaks = number_ticks(11))  +
    scale_y_continuous(limits = c(-1*ylimit, ylimit), breaks = number_ticks(11))  +
    geom_abline(intercept = ols$coefficients[1],
                slope = ols$coefficients[2])
}


volcanoplot_limma <- function(DEG_data_set, pcutoff=0.05, size=1){
  p <- ggplot(data= DEG_data_set, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
    geom_point(alpha=0.8, size = size) +
    theme_bw(base_size = 15) +
    theme(plot.title=element_text(hjust=0.5),   #  标题居中
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) + # 网格线设置为空白
    geom_hline(yintercept=-log10(pcutoff) ,linetype=4) +
    scale_color_manual(name = "", 
                       values = c("red", "green", "gray"),
                       limits = c("UP", "DOWN", "NOT"))
  return(p)
}


volcanoplot_deseq <- function(DEG_data_set){
  p <- ggplot(data= DEG_data_set, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
    geom_point(alpha=0.8, size = 2) +
    theme_bw(base_size = 15) +
    theme(plot.title=element_text(hjust=0.5),   #  标题居中
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) + # 网格线设置为空白
    geom_hline(yintercept=-log10(pcutoff) ,linetype=4) +
    #geom_vline(xintercept=c(-1,1) ,linetype=4 ) +
    scale_color_manual(name = "", 
                       values = c("red", "green", "gray"),
                       limits = c("UP", "DOWN", "NOT"))
  return(p)
}


pcaplotting <- function(pcaData, color=pcaData$condition, file="figures/pca.pdf"){
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=color)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
  ggsave(file,width = 10,height = 8)
}


# 读取猴子全脑数据的一个函数
read_monkey_data_and_region_rename <- function(monkeypath='D:/Monkey_Data/'){
  d99atlas        <- readNIfTI(paste0(monkeypath, 
                                     "Region_info_20220826/D99_2ion_combined173_final.nii.gz"))
  sample_metadata <- read.csv(paste0(monkeypath,
                                     "Sample_info_NCX_WholeBrain_GeneCounts_20220826/Information of sequenced samples_update_full878_filter819.CSV"), header=T)
  d99_metadata    <- read_excel(paste0(monkeypath,
                                     "Region_info_20220826/D99_2ion_combined173_CentroidCoordinate_0p25mm_20201112.xlsx"))
  # 处理脑残样本脑区命名
  wrong_region_name <- unique(sample_metadata$roi173)[unique(sample_metadata$roi173) %in% d99_metadata$Area==F]
  # 改样本表格里的脑残roi173命名
  sample_metadata[sample_metadata$roi173=="7m(PGm)",]$roi173 <- '7m'
  sample_metadata[sample_metadata$roi173=="F1_(4)",]$roi173 <- 'F1'
  sample_metadata[sample_metadata$roi173=="GPeGPi",]$roi173 <- 'GPe/GPi'
  sample_metadata[sample_metadata$roi173=="8Av",]$roi173 <- '8A'
  sample_metadata[sample_metadata$roi173=="7a_(Opt/PG)",]$roi173 <- '7a'
  sample_metadata[sample_metadata$roi173=="7b_(PFG/PF)",]$roi173 <- '7b'
  # 发现注释表格中RM94脑区有错
  ## 1. RM94的丘脑为Tha而不是TH，D99中的丘脑亚区也应标为Tha
  sample_metadata[sample_metadata$rm94=="TH",]$rm94 <- 'Tha'
  sample_metadata[sample_metadata$roi173=="CM",]$rm94 <- 'A2'
  ## 2. A2、PFCcl、CCr、HT没有采样
  # 还剩下STG和TG没有处理
  # 返回一个list
  return(list(d99atlas=d99atlas, sample_metadata=sample_metadata, d99_metadata=d99_metadata))
}


# 处理脑区对应value的一个函数
region_name2region_value <- function(region_name, hemi, d99_metadata){
  # 在剩下没处理的那两个脑残命名之外的就直接对应
  if (region_name != 'STG' & region_name != 'TG') {
    return(d99_metadata[d99_metadata$Area==region_name & d99_metadata$`Brain(L/R)`==hemi,]$Value)
  }
  else if (region_name=='STG' & hemi=='L') {
    return(c(101,68,15))
  }
  else if (region_name=='STG' & hemi=='R') {
    return(c(326,293,240))
  }
  else if (region_name=='TG' & hemi=='L') {
    return(c(59,154,78,110,135,140,115))
  }
  else if (region_name=='TG' & hemi=='R') {
    return(c(284,379,303,335,360,365,340))
  }
}


# 处理脑区对应value的一个函数，但是返回str
region_name2region_returnstr <- function(region_name, hemi, d99_metadata){
  # 在剩下没处理的那两个脑残命名之外的就直接对应
  if (region_name != 'STG' & region_name != 'TG') {
    return(as.character(d99_metadata[d99_metadata$Area==region_name & d99_metadata$`Brain(L/R)`==hemi,]$Value))
  }
  else if (region_name=='STG' & hemi=='L') {
    return('101,68,15')
  }
  else if (region_name=='STG' & hemi=='R') {
    return('326,293,240')
  }
  else if (region_name=='TG' & hemi=='L') {
    return('59,154,78,110,135,140,115')
  }
  else if (region_name=='TG' & hemi=='R') {
    return('284,379,303,335,360,365,340')
  }
}


# 样本对应到脑区value生成表格（非原始d99
sample2region_dataframe <- function(sample_metadata, d99_metadata){
  sample2region <- data.frame(sample=sample_metadata$No., region=sample_metadata$roi173, Side=sample_metadata$Side)
  value <- c()
  for (i in rownames(sample2region)){
    value[i] <- region_name2region_returnstr(sample2region[i,]$region, sample2region[i,]$Side, d99_metadata)
  }
  sample2region$value <- value
  return(sample2region)
}


# 生成d99脑区矩阵
generate_d99_exprSet <- function(exprSet,sample_metadata,d99_metadata){
  sample_metadata <- sample_metadata[sample_metadata$No. %in% colnames(exprSet),]
  sample2region   <- sample2region_dataframe(sample_metadata,d99_metadata)
  exprRegion <- data.frame(matrix(nrow = nrow(exprSet),ncol = 0))
  for (hemi in c("L","R")){
    for (each_region in unique(sample2region$region)){
      print(each_region)
      vectr <- rowMeans(select(exprSet,
                               sample2region[sample2region$region==each_region & sample2region$Side==hemi,]$sample))
      exprRegion[paste0(each_region,"_",hemi)] <- vectr
    }
  }
  rownames(exprRegion) <- rownames(exprSet)
  #删除含有缺失值的列
  library(dplyr)
  exprRegion <- exprRegion %>% select_if(~!any(is.na(.)))
  return(exprRegion)
}


# 读取筛选猴子皮层的一个函数
take_just_cortex_sample <- function(monkeypath='D:/Monkey_Data/'){
  vsd_cortex <- read.table(paste0(monkeypath, "Sample_info_NCX_WholeBrain_GeneCounts_20220826/mfas5_757samples_23613genes_vsd4_rmbatch.xls"), 
                           header=T, check.names = F)
  returnlist      = read_monkey_data_and_region_rename(monkeypath)
  sample_metadata = returnlist$sample_metadata
  d99_metadata    = returnlist$d99_metadata
  d99atlas        = returnlist$d99atlas
  #筛选皮层metadata
  sample_metadata <- sample_metadata[sample_metadata$No. %in% colnames(vsd_cortex),]
  return(list(sample_metadata=sample_metadata, d99_metadata=d99_metadata, vsd=vsd_cortex, d99atlas=d99atlas))
}

# 读取筛选猴子所有样本的一个函数
take_all_sample <- function(monkeypath='D:/Monkey_Data/'){
  vsd <- read.table(paste0(monkeypath, "Sample_info_NCX_WholeBrain_GeneCounts_20220826/mfas5_819samples_23605genes_vsd4_rmbatch.xls"), 
                           header=T, check.names = F)
  returnlist      = read_monkey_data_and_region_rename(monkeypath)
  sample_metadata = returnlist$sample_metadata
  sample_metadata = sample_metadata[sample_metadata$No. %in% colnames(vsd),]
  d99_metadata    = returnlist$d99_metadata
  d99atlas        = returnlist$d99atlas
  #筛选皮层metadata
  #sample_metadata <- sample_metadata[sample_metadata$No. %in% colnames(vsd_cortex),]
  return(list(sample_metadata=sample_metadata, d99_metadata=d99_metadata, vsd=vsd, d99atlas=d99atlas))
}


# 计算脑区和roi掩膜之间交集占脑区比例的一个函数 (傻方法)
region_mask_overlap_ratio <- function(region_value, roi_mask, d99atlas){
  region <- niftiarr(d99atlas, d99atlas %in% region_value)
  overlap <- niftiarr(region, region*roi_mask)
  overlap_sum <- sum(overlap, na.rm=T)
  region_sum <- sum(region, na.rm=T)
  result <- overlap_sum/region_sum
  return(result)
}


atlas_table <- function(atlas){
  t <- table(atlas)
  t <- as.data.frame(t)
  t <- t[-1,]
  colnames(t) <- c('atlas', 'Freq')
  return(t)
}


# 聪明办法
calculate_d99_roi_cross_ratio_in_one_hemi <- function(sample_metadata, d99_metadata, d99atlas, d99table, roi, hemi, raw=F, plot=T){
  if(raw){
    cross <- niftiarr(d99atlas, roi*d99atlas)
    crosstable <- atlas_table(cross)
    ratio_in_roi = c()
    for (region in unique(sample_metadata$roi173)){
      region_value = region_name2Raw_region_value(region, hemi, d99_metadata)
      ratio        <- sum(crosstable[crosstable$atlas %in% region_value, 'Freq'])/sum(d99table[d99table$atlas %in% region_value,'Freq'])
      ratio_in_roi <- append(ratio_in_roi, ratio)
    }
  }
  else{
    cross <- niftiarr(d99atlas, roi*d99atlas)
    crosstable <- atlas_table(cross)
    ratio_in_roi = c()
    for (region in unique(sample_metadata$roi173)){
      region_value = region_name2region_value(region, hemi, d99_metadata)
      ratio        <- sum(crosstable[crosstable$atlas %in% region_value, 'Freq'])/sum(d99table[d99table$atlas %in% region_value,'Freq'])
      ratio_in_roi <- append(ratio_in_roi, ratio)
    }
  }
  if(plot){
    plot(sort(ratio_in_roi))
  }
  return(ratio_in_roi)
}


# 分别对半脑的每一个脑区循环计算脑区roi交集比例的一个函数 （傻方法）
calculate_overlap_ratio_for_each_region_in_one_hemi <- function(sample_metadata, d99_metadata, d99atlas, roi, hemi, raw=F){
  if (raw){
    ratio_in_roi = c()
    for (region in unique(sample_metadata$roi173)){
      region_value = region_name2Raw_region_value(region, hemi, d99_metadata)
      ratio_in_roi <- append(ratio_in_roi, region_mask_overlap_ratio(region_value, roi, d99atlas))
    }
  }
  else{
    ratio_in_roi = c()
    for (region in unique(sample_metadata$roi173)){
      region_value = region_name2region_value(region, hemi, d99_metadata)
      ratio_in_roi <- append(ratio_in_roi, region_mask_overlap_ratio(region_value, roi, d99atlas))
    }
  }
  plot(sort(ratio_in_roi))
  return(ratio_in_roi)
}


# 样本对应到原始d99脑区value的函数，但是返回字符串
region_name2Raw_region_str <- function(region_name, hemi, d99_metadata){
  # 在剩下没处理的那两个脑残命名之外的就直接对应
  if (region_name != 'STG' & region_name != 'TG' & region_name != 'amy' 
      & region_name != 'amy' & region_name != '10m' & region_name != 'HC' & region_name != 'V2' & region_name != '5' & region_name != '8A' & region_name != '45' & region_name != 'LGN' & region_name != 'MD') {
    return(as.character(d99_metadata[d99_metadata$Area==region_name & d99_metadata$`Brain(L/R)`==hemi,]$Value))
  }
  else if (region_name=='STG' & hemi=='L') {
    return('101,68,15')
  }
  else if (region_name=='STG' & hemi=='R') {
    return('326,293,240')
  }
  else if (region_name=='TG' & hemi=='L') {
    return('59,154,78,110,135,140,115')
  }
  else if (region_name=='TG' & hemi=='R') {
    return('284,379,303,335,360,365,340')
  }
  else if (region_name=='amy' & hemi=='L') {
    return('4,13,21,33,42,58,138,147,168,179')
  }
  else if (region_name=='amy' & hemi=='R') {
    return('229,238,246,258,267,283,363,372,393,404')
  }
  else if (region_name=='10m' & hemi=='L') {
    return('14,47')
  }
  else if (region_name=='10m' & hemi=='R') {
    return('239,272')
  }
  else if (region_name=='HC' & hemi=='L') {
    return('192,189,157,54,83,19,88,190,191')
  }
  else if (region_name=='HC' & hemi=='R') {
    return('417,414,382,279,308,244,313,415,416')
  }
  else if (region_name=='V2' & hemi=='L') {
    return('131,172,174')
  }
  else if (region_name=='V2' & hemi=='R') {
    return('356,397,399')
  }
  else if (region_name=='5' & hemi=='L') {
    return('43,77,134')
  }
  else if (region_name=='5' & hemi=='R') {
    return('268,302,359')
  }
  else if (region_name=='8A' & hemi=='L') {
    return('148,51')
  }
  else if (region_name=='8A' & hemi=='R') {
    return('376,276')
  }
  else if (region_name=='45' & hemi=='L') {
    return('124,143')
  }
  else if (region_name=='45' & hemi=='R') {
    return('349,368')
  }
  else if (region_name=='LGN' & hemi=='L') {
    return('200,201')
  }
  else if (region_name=='LGN' & hemi=='R') {
    return('425,426')
  }
  else if (region_name=='MD' & hemi=='L') {
    return('206,209')
  }
  else if (region_name=='MD' & hemi=='R') {
    return('431,434')
  }
}


# 样本对应到原始d99脑区value的函数
region_name2Raw_region_value <- function(...){
  return(as.integer(strsplit(region_name2Raw_region_str(...),split = ',')[[1]]))
}


# 样本对应到脑区value生成表格（原始d99
sample2_Raw_region_dataframe <- function(sample_metadata, d99_metadata){
  sample2region <- data.frame(sample=sample_metadata$No., region=sample_metadata$roi173, Side=sample_metadata$Side)
  value <- c()
  for (i in rownames(sample2region)){
    value[i] <- region_name2Raw_region_str(sample2region[i,]$region, sample2region[i,]$Side, d99_metadata)
  }
  sample2region$value <- value
  return(sample2region)
}

process_time <- function(fun_name){
  t1=proc.time()
  a <- fun_name
  t2=proc.time()
  t=t2-t1
  print(paste0('processing time：',t[3][[1]],'秒'))
  return(a)
}


# 将nii脑区分区图像右脑加上最大值的函数
nii_right_plus <- function(nii,plus){
  mask <- niftiarr(nii, nii>0)
  nii@.Data[0:ceiling(dim(nii)[1]/2),,] = nii@.Data[0:ceiling(dim(nii)[1]/2),,] + plus
  nii <- niftiarr(nii, nii*mask)
  return(nii)
}


# 获得kegg通路的所有基因名(需要包KEGGREST调用kegg的api)
get_kegg_genesymbol <- function(keggid){
  gs <- keggGet(keggid)
  genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
  genesymbol<- genes[1:length(genes)%%3 ==2]
  print(gs[[1]]$NAME)
  return(genesymbol)
}


correlationplot <- function(x, y, xlabel='x', ylabel='y'){
  corrdata <- data.frame(x=x, y=y)
  ggplot(corrdata, aes(x, y)) + 
    xlab(xlabel)+ylab(ylabel)+
    geom_point(size=2.5)+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'pearson', aes(x =x, y =y),  label.x.npc = "middle",label.y.npc = "bottom",)+ 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}