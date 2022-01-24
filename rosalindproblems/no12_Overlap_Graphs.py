'''
Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.
Return: The adjacency list corresponding to O3. You may return edges in any order.
'''

from Bio import SeqIO
path = r'C:\Users\11262\PycharmProjects\rosalindproblems\no12.fasta'
for seqrecord1 in SeqIO.parse(path,'fasta'):
    for seqrecord2 in SeqIO.parse(path,'fasta'):
        if seqrecord1.name != seqrecord2.name and seqrecord1.seq[-3:] == seqrecord2.seq[:3]:
            print(seqrecord1.name,seqrecord2.name)