MonkGroupData <- "D:/Myworks/research_in_monkeyhuman_transcriptome/ImagingData/Monk_S63_GroupAvg"
HCPGroupData <- "D:/Myworks/research_in_monkeyhuman_transcriptome/ImagingData/HCP_S1200_GroupAvg_v1"


# 皮层xifti数据处理相关
# 皮层分割函数
# trim 截断平均
parcAverageSurfaceData <- function(xii_Data, xii_Parc, trim=0) {
  # 注意分割向量的长度和皮层数据向量长度是否一致，如果不一致。
  # 在使用该函数之前，需要remove_from_mwall
  # 直接使用parcel数据比使用parc“xifti”对象更容易。因此，从将分割转化为向量开始
  parc_vec <- c(as.matrix(xii_Parc))
  # 对于每个parcel，将获得其组成顶点的平均值------------------------
  Data_mat  <- as.matrix(xii_Data)
  n_parc <- tail(xii_Parc$meta$cifti$labels[[1]]$Key, 1)
  Data_pmean_mat <- matrix(nrow=n_parc, ncol=ncol(Data_mat))
  for (p in 1:n_parc) {
    data_p <- Data_mat[parc_vec==p,]
    Data_pmean_mat[p, ] <- apply(data_p, 2, mean, trim=trim, na.rm=T)
  } 
  xii_Data_pmean <- rbind(NA, Data_pmean_mat)[parc_vec + 1,]
  print(dim(xii_Data_pmean))
  print(dim(Data_pmean_mat))
  xii_pmean <- newdata_xifti(xii_Data, xii_Data_pmean)
  return(list(Parc_xii=xii_pmean, Parc_matrix=Data_pmean_mat))
}


# 根据分区矩阵生成全脑xii
# 注意；这里是根据xii_Parc替换生成新的cifiti文件，所以列不同，
# 所以会损失列名的信息，在meta属性中，主要是sample names，比如100206_thickness之类的
# 但sample的顺序是一定的，所以还是可以根据顺序的索引获得对应被试的信息
# df needs to be RxS (regions x samples)
parcMean2xiiData <- function(df, xii_Parc) {
  parc_vec <- c(as.matrix(xii_Parc))
  n_parc <- tail(xii_Parc$meta$cifti$labels[[1]]$Key, 1)
  xii_Data_pmean <- rbind(NA, df)[parc_vec + 1,]
  xii_pscalar <- newdata_xifti(ciftiTools:::convert_to_dscalar(xii_Parc), xii_Data_pmean)
  return(xii_pscalar)
}




# 画图相关
# 相关性分析及画图
correlationplot <- function(x, y, xlabel='x', ylabel='y', title="", file="", width=8, height=8){
  corrdata <- data.frame(x=x, y=y)
  ggplot(corrdata, aes(x, y)) + 
    xlab(xlabel)+ylab(ylabel)+
    geom_point(size=2.5)+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'pearson', aes(x=x, y=y),  label.x.npc = "middle",label.y.npc = "bottom",)+ 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle(title) -> p
  if (file != "") {
    ggsave(p, file=file, width=width, height=height)
  }
  return(p)
}

# 保存热图
save_pheatmap_pdf <- function(x, filename, width=12, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}













# others
