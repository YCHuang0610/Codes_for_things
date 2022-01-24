from Bio import SeqIO
from Bio.SeqUtils import GC
path = r'C:\Users\11262\PycharmProjects\rosalindproblems\no6.fasta'
gc=[]
seqrecords = SeqIO.parse(path,'fasta')
for seqrecord in seqrecords:
    gc.append(GC(seqrecord.seq))
max = max(gc)
i = gc.index(max)
list = list(SeqIO.parse(path,'fasta'))
print(list[i].id,max)


