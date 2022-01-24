from Bio import SeqIO
record = SeqIO.read('no21.fasta','fasta')
seq = record.seq
'''
找回文序列
制酶识别DNA序列中的回文序列
挨个扫描序列中长度在4-12之间的片段，与反向互补序列比较，相同即为酶切位点序列
'''
seq_com = seq.complement()
print(seq)
print(seq_com)
for i in range(0,len(seq)):
    for j in range(4,13):
        if (i+j)>len(seq):
            break
        else:
            a = seq[i:i+j]
            if a[::-1] == seq_com[i:i+j]:
                print(i+1,j)