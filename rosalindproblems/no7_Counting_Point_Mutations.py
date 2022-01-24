from Bio.Seq import Seq
s = Seq('GAGCCTACTAACGGGAT')
t = Seq('CATCGTAATGACGGCCT')
hamming = 0
for i in range(len(s)):
    if s[i] != t[i]:
        hamming += 1
print(hamming)

