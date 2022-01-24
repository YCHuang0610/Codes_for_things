from Bio.Seq import Seq
s = Seq('GATATATGCATATACTT')
t = Seq('ATAT')
for i in range(len(s)-len(t)):
    if s[i:i+len(t)] == t:
        print(i+1)
#也可以用re模块