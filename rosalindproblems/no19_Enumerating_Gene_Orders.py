'''
全排列 比如输入123
输出123 132 213 231 312 321
'''
a = (input('输入：')).split(',')
#全排列算法
#下降二叉树法
def perm(seq,start,end):
    global num
    if start == end:  #当前n-1个数都完成交换之后输出
        print(''.join(seq))
        num += 1
    else:
        for i in range(start,end):#将第一个数与每个数进行交换
            seq[i],seq[start] = seq[start],seq[i]
            perm(seq,start+1,end)
            seq[i], seq[start] = seq[start], seq[i]
num = 0
perm(a,0,len(a))
print(num)