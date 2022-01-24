from Bio import SeqIO
rec_list = list(SeqIO.parse('no22.fasta','fasta'))
seq = rec_list[0].seq
intron1 = rec_list[1].seq
intron2 = rec_list[2].seq
index1 = seq.index(intron1)
index2 = seq.index(intron2)
print(index1,index2)
newseq = seq[0:index1]+seq[index1+len(intron1):index2]+seq[index2+len(intron2)::]
pro = newseq.translate()
print(pro)

