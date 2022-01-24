from sklearn import datasets
import matplotlib.pyplot as plt
import numpy as np


# 二分类问题，三层神经元，神经元个数分别是2,4,1

# 添加偏置神经元的函数
def plus_bias(X):
    x_0 = np.full((X.shape[0], 1), 1)
    return np.column_stack((X, x_0))

# 激活函数
def sigmod(x):
    return 1/(1+np.exp(-x))

# 激活函数的导数
def sigmod_derivative(x):
    return x*(1-x)



class NN():
    def __init__(self, X, y):
        self.input = X  #200*2
        self.y = y #200*1
        self.m = X.shape[0]
        self.inputlayer = X.shape[1]
        self.hiddenlayer = 4
        self.outputlayer = 1
        self.w_ih = np.random.rand(self.inputlayer + 1, self.hiddenlayer) #3*4
        self.w_ho = np.random.rand(self.hiddenlayer + 1, self.outputlayer)  #5*1
        
    # 前向传播
    def forword(self):
        self.a_hidden = sigmod(np.dot(plus_bias(self.input), self.w_ih)) #200,4
        self.output = sigmod(np.dot(plus_bias(self.a_hidden),self.w_ho)) #200*1
        
    # 反向传播，求导数
    def back(self):
        # 输出层误差 δ
        out_error = self.output - self.y #200*1
        # 链式法则计算中间层与输出层间权重的导数
        # 先将误差传播到激活前：用误差乘以激活函数的导数 （越接近的误差越小）
        # 用中间层神经元的值的转置乘传到激活前的误差 （输入值越大权重越大）
        self.d_w_ho = np.dot(plus_bias(self.a_hidden).T, out_error*sigmod_derivative(self.output)) #5*1
        # 用误差乘权重，将误差传到中间层
        hidden_error = np.dot(out_error, self.w_ho[0:4].T) #200,4
        # 同上，用链式法则求输入层与中间层权重的导数
        self.d_w_ih = np.dot(plus_bias(self.input).T, hidden_error*sigmod_derivative(self.a_hidden))
        
    # 训练,学习率0.001,迭代一百万次
    def train(self, lr = 0.001, epochs = 1000000):
        for epoch in range(epochs):
            # 先前向传播再反向传播
            self.forword()
            self.back()
            # 梯度下降更新层间权重
            self.w_ih -= lr*self.d_w_ih
            self.w_ho -= lr*self.d_w_ho
    # 用学习的权重预测函数
    def predict(self, X):
        a_hidden = sigmod(np.dot(plus_bias(X), self.w_ih))
        output = sigmod(np.dot(plus_bias(a_hidden),self.w_ho))
        print(output)
        return np.where(output >= 0.5, 1, 0)



# 抄的一段画决策边界的函数
def plot_decision_boundary(pred_func):
    # 设定最大最小值，附加一点点边缘填充
    x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
    y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5
    h = 0.01

    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

    # 用预测函数预测一下
    Z = pred_func(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)

    # 然后画出图
    plt.contourf(xx, yy, Z, cmap=plt.cm.Spectral)
    plt.scatter(X[:, 0], X[:, 1], c=y, cmap=plt.cm.Spectral)




if __name__ == '__main__':
    np.random.seed(0)
    
    # 月牙形加噪音创造散点数据
    X,y = datasets.make_moons(200,noise=0.2)
    y_b = y.reshape((200,1))
    # 画散点图
    plt.scatter(X[:,0],X[:,1],c=y)
    
    #训练神经网络
    neural = NN(X,y_b)
    neural.train()

    # 画决策边界
    plot_decision_boundary(lambda x: neural.predict(x))
    plt.show()
