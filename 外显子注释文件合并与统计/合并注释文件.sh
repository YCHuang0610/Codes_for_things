#!/bin/sh 
# TR硬盘目录
wk=".....Variation/" #填入工作目录
cd "$wk"
usr=~/Desktop/hyc/TR_joint/ #填入输出目录
# find all annotated snp
find . -name '*.samtools.snp.annovar.hg19_*' -print > $usr/TR.snp.config.txt

#snp
# joint snp select H priority
for i in $(cat $usr/TR.snp.config.txt);do
	sample=$(echo $i | awk -F '/' '{print $2}')
	zcat < $i | grep ^H | sed "s/$/\t${sample}/" >> $usr/TR.snp.H.joint.tsv
done

# 将列名放到header.txt里面存着
# zcat < "$wk"/TAR_760/SNP/TAR_760.samtools.snp.annovar.hg19_multianno.xls.gz | head -n 1 > $usr/header.txt
# 手动添加一个列名在最后叫SampleID
header=$(cat $usr/header.txt)
sed -i "" "1s/^/${header}\n/" $usr/TR.snp.H.joint.tsv

#indel
find . -name '*.samtools.indel.annovar.hg19_*' -print > $usr/TR.indel.config.txt
for i in $(cat $usr/TR.indel.config.txt);do
	sample=$(echo $i | awk -F '/' '{print $2}')
	zcat < $i | grep ^H | sed "s/$/\t${sample}/" >> $usr/TR.indel.H.joint.tsv
done

header=$(cat $usr/header.txt)
sed -i "" "1s/^/${header}\n/" $usr/TR.indel.H.joint.tsv
