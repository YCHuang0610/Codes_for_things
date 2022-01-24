from Bio import SeqIO
def align_count(seq1,seq2):     #首尾重叠数目
    list_i =[]
    for i in range(0,len(seq1)):
        if seq2[0:i] == seq1[-i::]:
            list_i.append(i)
    try:
        return max(list_i)
    except:
        return 0
def orders(seq1,seq2):      #我们规定前面的尾对着后面的头，
    flag1 = align_count(seq1,seq2)
    flag2 = align_count(seq2,seq1)
    if flag1>=flag2:
        return flag1,1 #1在前2在后
    else:
        return flag2,0
def addseq(n,seq1,seq2):     #前面的尾巴和后面的头相同，将两个序列重叠合并
    newseq = seq1[:-n]+seq2
    return newseq

if __name__ == '__main__':
    recordlist = list(SeqIO.parse('no25.fasta', 'fasta'))
    seqlist = [record.seq for record in recordlist]
    while len(seqlist)>1:      #两两比对后取重叠部分最多那俩，然后从原有的序列列表里将这两个序列变成他们的重叠序列，然后剩下的再来一遍，一直到只剩下一个序列为止
        n_list = []   #储存每次的重叠数目
        newseqlist = [] #储存每次的合并重叠序列
        i_list = []
        j_list = []
        #这四个列表对每组数据的索引是一样的
        for i in range(0,len(seqlist)-1):
            for j in range(i+1,len(seqlist)):
                n,front = orders(seqlist[i],seqlist[j])
                n_list.append(n)
                i_list.append(i)
                j_list.append(j)
                if front == 1:
                    newseq = addseq(n,seqlist[i],seqlist[j])
                    newseqlist.append(newseq)
                if front == 0:
                    newseq = addseq(n,seqlist[j],seqlist[i])
                    newseqlist.append(newseq)
        maxflag = max(n_list)   #取重叠部分最多的
        maxindex = n_list.index(maxflag)
        i = i_list[maxindex]
        j = j_list[maxindex]
        seqlist[i] = newseqlist[maxindex]
        seqlist.pop(j)
    print(seqlist[0])


