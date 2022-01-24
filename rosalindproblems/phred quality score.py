'''
质量值是测序错误率的对数（10为底数）乘以-10（并取整）
描述的是每个测序碱基的可靠程度，用ASCII码表示。
'''
qual = 'JJJJJIIJJJJJJHIHHHGHFFFFFFCEEEEEDBD'
phred_score = [ord(q)-33 for q in qual] #33是测序平台质量体系所加的值Phred33
print(phred_score)
mistake = [10**(-q/10.0) for q in phred_score]
print(mistake)
