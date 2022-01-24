from Bio.Seq import Seq
s = Seq('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
p = s.translate()
print(p)