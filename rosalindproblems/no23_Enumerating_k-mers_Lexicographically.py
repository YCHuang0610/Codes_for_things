seq = 'ACGT'
n = 0
for i in range(0,len(seq)):
    for j in range(0,len(seq)):
        print(seq[i]+seq[j])
        n += 1
print(n)

