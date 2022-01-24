import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
path = r'C:\Users\11262\PycharmProjects\rosalindproblems\no10.fasta'
seq_array = np.array([list(seqrecord) for seqrecord in SeqIO.parse(path,'fasta')])
print(seq_array)
print(seq_array.shape)
a = np.sum(seq_array=='A',axis=0) #求每一列的A的数目
print('A:',a)
t = np.sum(seq_array=='T',axis=0)
print('T:',t)
g = np.sum(seq_array=='G',axis=0)
print('G:',g)
c = np.sum(seq_array=='C',axis=0)
print('C:',c)
array = np.array([a,t,g,c])
newseqindex = np.argmax(array,axis=0) #np.argmax返回最大值对应索引，储存索引列表，0是a，1是t，2是g，3是c
replacediclist = 'ATGC'               #替换字典列表其索引正好对应着0是a，1是t，2是g，3是c
newseq = []                           #储存新序列的列表
for i in range(len(newseqindex)):
    newseq.append(replacediclist[newseqindex[i]])
seq = ''.join(newseq)                 #将列表转换为字符串
print(seq)


