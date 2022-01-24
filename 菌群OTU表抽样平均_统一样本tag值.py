'''
OTU表抽平
'''
import pandas as pd
import numpy as np
from collections import Counter


df = pd.read_table(r'C:\Users\11262\Downloads\EasyAmplicon-master\EasyAmplicon-master\otutab.txt',header=0,index_col=0)
#设置最大tag数，小于原表格最小值。
size = 32086
print(df)

def flattt(df):
    #将每一列按照reads次数扁平化为一维
    dict = {}
    for col in df.columns.tolist():
        s = df[col]
        list = np.repeat(s.index.tolist(),s.values)
        dict[col] = list
    return dict

dict = flattt(df)

def choice_count(list,size):
    #对每列不放回抽样，并统计抽样后每个元素出现次数
    list_choice = np.random.choice(list,size,replace=False, p=None)
    dict_count = Counter(list_choice)
    return dict_count

#将抽样结果变回dataframe
df_dict = {}
for col in df.columns.tolist():
    df_dict[col] = pd.Series(choice_count(dict[col],size),index=df._stat_axis) #dataframe的行名是df._stat_axis
df_rare = pd.DataFrame(df_dict)
df_rare = df_rare.fillna(0)
df_rare = df_rare.astype(int)

print(df_rare)

df_rare.to_csv(r'C:\Users\11262\Downloads\EasyAmplicon-master\EasyAmplicon-master\otutab_rare.txt', sep='\t')
