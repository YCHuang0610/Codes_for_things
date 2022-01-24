import numpy as np

# 假设方程
def hypoth(Xb,thetas):
    return Xb.dot(thetas)


class LinerRegression(object):
    def __init__(self):
        pass
# 拟合
    def fit(self,X,y,alpha=0.01,accuracy=0.0001):
        self.m = X.shape[0]
        self.n = X.shape[1]
        dimension = self.n + 1
        # 初始化theta
        self.thetas = np.full(dimension, 0.5)
        # 添加偏置
        x_0 = np.full((self.m, 1), 1)
        Xb = np.column_stack((x_0, X))
        # 梯度下降
        count = 1
        while True:
            oldJ = self.costFunc(Xb, y)
            # 同步更新thetas
            temp = np.zeros(self.thetas.shape) #用temp存储中间变量，实现同步更新
            for j in range(dimension):
                sum = 0
                for i in range(self.m):
                    sum += (hypoth(Xb[i,:],self.thetas) - y[i])*Xb[i,j] #求导数梯度
                temp[j] = self.thetas[j] - alpha/self.m*sum #下降
            self.thetas = temp #更新thetas
            newJ = self.costFunc(Xb, y)
            # 梯度下降结束
            if newJ == oldJ or abs(newJ - oldJ) <= accuracy:
                print('*****梯度下降完了*****')
                print('收敛到代价函数值为', newJ)
                print('迭代次数', count)
                print('***********')
                break
            print('迭代第%d次' % count)
            count += 1
            
    # 损失函数
    def costFunc(self,Xb,y):
        cost = 0.0
        for i in range(self.m):
            cost += (hypoth(Xb[i,:], self.thetas) - y[i])**2
        return cost/2/self.m
    
    # 预测值
    def predict(self,X):
        x_0 = np.full((X.shape[0], 1), 1)
        Xb = np.column_stack((x_0, X))
        return hypoth(Xb, self.thetas)

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    X = 2 * np.random.rand(100, 1)
    # y=3x+4 数据
    y = 4 + 3 * X + np.random.randn(100, 1)
    plt.scatter(X,y)
    # 初始化线性回归类
    Line = LinerRegression()
    # 拟合
    Line.fit(X,y)
    # 画拟合直线
    X_new = np.array([[0], [2]])
    y_predict = Line.predict(X_new)
    plt.plot(X_new,y_predict,color='red')
    plt.show()
