'''
Given: A collection of k (k<100) DNA strings of length at most 1 kbp, each in FASTA format.
Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)
'''
path = r'C:\Users\11262\PycharmProjects\rosalindproblems\no14.fasta'
from Bio import SeqIO
import re
def motif(seq):    #找子序列的函数，返回一个子序列列表
    list = []
    for i in range(len(seq) - 1):  # 子序列从第几位开始
        for j in range(2,len(seq)+1-i):  #子序列长度，大于2有意义
            list.append(seq[i:i+j])
    return list
def shortest(seqpath):   #返回最短的序列和序列的索引
    lens = []
    for seqrecord in SeqIO.parse(seqpath,'fasta'):
        lens.append(len(seqrecord))
    i = lens.index(min(lens))
    seqlist = list(SeqIO.parse(seqpath,'fasta'))
    return seqlist[i].seq,i

if __name__ == '__main__':
    seqlist = list(SeqIO.parse(path,'fasta'))
    seq, a = shortest(path)
    for s in motif(seq):
        flag = 0
        for i in range(len(seqlist)):
            if i == a:
                flag += 1
                if flag == len(seqlist):
                    print(s)
            else:
                if re.search(str(s),str(seqlist[i].seq)):
                    flag += 1
                    if flag == len(seqlist):
                        print(s)







