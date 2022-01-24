'''
不知道哪个是
'''

from Bio import SeqIO

path = r'C:\Users\11262\PycharmProjects\rosalindproblems\no18.fasta'
record = SeqIO.read(path,"fasta")
a = record.seq
b = a.reverse_complement()

dna = [a,a[1:],a[2:],b,b[1:],b[2:]]

for seq in dna:
    try:
        pro = seq.translate(to_stop =True)
        ilist = [i for i,x in enumerate(pro) if x=='M']
        for i in ilist:
            print(pro[i:])
    except:
        pass



