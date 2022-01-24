import numpy as np
import math
from sklearn import datasets
from collections import Counter


infinity = float(-2**31)

def sigmodformatrix(Xb, thetas):
    t = - Xb.dot(thetas)
    s = np.zeros(t.shape[0])
    for i in range(len(s)):
        s[i] = 1 / (1 + math.exp(t[i]))
    return s
  
def sigmod(Xi,thetas):
    t = -np.sum(Xi*thetas) #-Xi.dot(thetas)
    s = 1/(1 + math.exp(t))
    return s

class Logisticregression():
    thetas = None
    m = 0
    def __init__(self):
        pass

    def fit(self,X,y,alpha=0.01,accuracy=0.00001):
        self.n = X.shape[1]
        self.m = X.shape[0]
        dimension = self.n + 1
        self.thetas = np.full(dimension,0.5)
        x_0 = np.full((self.m,1),1)
        # 添加x_0=1的一列
        Xb = np.column_stack((x_0,X))
        # 梯度下降
        count = 1
        while True:
            oldJ_theta = self.costFunc(Xb,y)
            #同步更新theta
            temp = np.zeros(self.thetas.shape)
            for j in range(dimension):
                temp[j] = self.thetas[j] - alpha*np.sum((sigmodformatrix(Xb,self.thetas) - y)*Xb[:,j])
                self.thetas = temp
            newJ_theta = self.costFunc(Xb,y)
            if newJ_theta == oldJ_theta or abs(newJ_theta-oldJ_theta) < accuracy:
                print('*****梯度下降完了*****')
                print('收敛到代价函数值为',newJ_theta)
                print('迭代次数',count)
                print('***********')
                break
            print('迭代第%d次' % count)
            count += 1
    def costFunc(self,Xb,y):
        cost = 0.0
        for i in range(self.m):
            h_theta = sigmod(Xb[i,],self.thetas)
            if h_theta == 1 or h_theta == 0:
                return infinity
            cost += y[i]*math.log(h_theta)+(1-y[i])*math.log(1-h_theta)
        return -1/self.m*cost
    def predict(self,X):
        x_0 = np.full((len(X),1),1)
        Xb = np.column_stack((x_0,X))
        s = sigmodformatrix(Xb,self.thetas)
        s[s>=0.5] = 1
        s[s<0.5] = 0
        return s.astype(np.int)

if __name__ == '__main__':
    from sklearn.model_selection import train_test_split

    iris = datasets.load_iris()
    X = iris['data']
    y = iris['target']
    X = X[y!=2]
    y = y[y!=2]
    X_train, X_text, y_train, y_text = train_test_split(X,y)
    Logistic = Logisticregression()
    Logistic.fit(X_train, y_train)
    y_predict = Logistic.predict(X_text)
    print('参数',Logistic.thetas)
    print(y_predict)
    print(y_text)
