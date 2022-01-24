#合并all统计各方面成表
snp <- read.csv("TR.snp.H.joint.tsv",header=T,sep='\t',na.string='.')
indel <- read.csv("TR.indel.H.joint.tsv",header=T,sep='\t',na.string='.')
snp$ExonicFunc[is.na(snp$ExonicFunc)] <- 'splicing SNV'
indel$ExonicFunc[is.na(indel$ExonicFunc)] <- 'splicing INDEL'
all <- rbind(snp,indel)

snp$RefAlt <- paste(snp$REF,snp$ALT,sep='>')
refalt <- with(snp,table(RefAlt))
refalt <- as.data.frame(refalt)
func <- with(all,table(Func))
func <- as.data.frame(func)
exonicfunc <- with(all,table(ExonicFunc))
exonicfunc <- as.data.frame(exonicfunc)
gene <- with(all,table(GeneName))
gene <- as.data.frame(gene)

#wordcloud
library(wordcloud)
wordcloud(gene$GeneName,gene$Freq,col=brewer.pal(8,"Dark2"),max.words = 50)

#变异类型统计
library(ggplot2)
p <- ggplot(data=exonicfunc,aes(x=ExonicFunc,y=Freq,
                                fill=ExonicFunc))+geom_bar(stat='identity'
                                                          )+coord_flip()+guides(fill=F
                                                                               )+theme_bw()
p
ggsave("exonicfunc.pdf",width = 6,height = 4)

#SNV类型
refalt <- refalt[order(refalt$Freq),]
#排序
refalt$RefAlt <- factor(refalt$RefAlt,levels = refalt$RefAlt)
#添加转换颠换类型
refalt$class <- c(rep("transversions",8),rep("transitions",4))
p <- ggplot(data=refalt,aes(x=RefAlt,y=Freq,
                             fill=class))+geom_bar(stat='identity'
                             )+coord_flip()+labs(x='SNV class',
                                                 y='Freq'
                                                )+theme_bw()+guides(fill=guide_legend(title = NULL)
                                                       )+theme(legend.justification=c(1,0), 
                                                               legend.position=c(0.9,0.1))
p
ggsave("SNVclass.pdf",width = 6,height = 3)

#TOP10GENE
gene <- gene[order(gene$Freq,decreasing = T),]
gene10 <- head(gene,10)
gene10$GeneName <- factor(gene10$GeneName,levels=gene10$GeneName)
p <- ggplot(data=gene10,aes(x=GeneName,
                            y=Freq,fill=GeneName))+geom_bar(stat='identity'
                                                           )+coord_flip()+guides(fill=F)+theme_bw()
p
ggsave("TOP10gene.pdf",width = 6,height = 3)

#瀑布图
library(GenVisR)
waterf <- data.frame(sample=all$SampleID,gene=all$GeneName,variant_class=all$ExonicFunc)
waterfall(waterf, fileType="Custom",maxGenes = 20,mainGrid=F,
          variant_class_order=c('missense SNV','frameshift deletion','splicing SNV','splicing INDEL',
                                'frameshift insertion','stopgain','stoploss',
                               'nonframeshift deletion','nonframeshift insertion',
                               'unknown') )
                           
