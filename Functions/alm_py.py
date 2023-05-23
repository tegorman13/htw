
# import sys
# !{sys.executable} -m pip install numpy

import numpy as np

def input_activation(x_target, c):
    return np.exp((-1 * c) * (x_target - inputNodes) ** 2)


def output_activation(x_target, weights, c):
    return np.dot(weights, input_activation(x_target, c))


def mean_prediction(x_target, weights, c):
    probability = output_activation(x_target, weights, c) / sum(output_activation(x_target, weights, c))
    return np.dot(outputNodes, probability)  # integer prediction


def update_weights(x_new, y_new, weights, c, lr, noise_sd=None):
    y_feedback_activation = np.exp(-1 * c * (y_new - outputNodes) ** 2)
    x_feedback_activation = output_activation(x_new, weights, c)
    return weights + lr * (y_feedback_activation - x_feedback_activation) * input_activation(x_new, c).T

def train_alm(dat, c=0.05, lr=0.5, weights=None):
    alm_train = np.empty(dat.shape[0])
    if weights is None:
        weights = np.empty(dat.shape[1])
        weights.fill(0.5)
    for i in range(dat.shape[0]):
        weights = update_weights(dat.input[i], dat.vx[i], weights, c, lr)
        resp = mean_prediction(dat.input[i], weights, c)
        alm_train[i] = resp
        weights[weights < 0] = 0
    return alm_train

def exam_prediction(x_target, weights, c,trainVec):
    #trainVec = sort(unique(x.learning))
    nearestTrain = trainVec[np.argmin(np.abs(trainVec-x_target))]
    aresp = mean_prediction(nearestTrain, weights, c)
    xUnder = nearestTrain if min(trainVec) == nearestTrain else trainVec[np.where(trainVec == nearestTrain) - 1]
    xOver = nearestTrain if max(trainVec) == nearestTrain else trainVec[np.where(trainVec == nearestTrain) + 1]
    mUnder = mean_prediction(xUnder, weights, c)
    mOver = mean_prediction(xOver, weights, c)
    exam_output = round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (x_target - nearestTrain), 3)
    return exam_output


def trainTest(dat, c=0.05, lr=0.5, weights, testVec, update_func, noise_sd):
    update_func = update_func
    alm_train = np.repeat(np.nan, dat.shape[0])
    for i in range(dat.shape[0]):
        weights = update_func(dat.input[i], dat.vx[i], weights, c, lr)
        resp = mean_prediction(dat.input[i], weights, c)
        alm_train[i] = resp
        weights[weights<0] = 0
    alm_pred = np.apply_along_axis(mean_prediction, 0, testVec, weights, c)
    exam_pred = np.apply_along_axis(exam_prediction, 0, testVec, weights, c, trainVec=np.array([1, np.sort(np.unique(dat.input))]))
    return {'almTrain': alm_train, 'almPred': alm_pred, 'examPred': exam_pred}



def gen_train(trainVec=[5,6,7],trainRep=3,noise=0):
    bandVec=[0,100,350,600,800,1000,1200]
    if(isinstance(trainVec,list)):trainVec=unlist(trainVec)
    ts = np.repeat(np.arange(1,len(trainVec)+1),trainRep)
    # print(trainVec)
    # print(len(ts)); print(len(trainRep))
    noiseVec=np.random.normal(0,1,len(ts))*noise
    if(noise==0) :noiseVec=noiseVec*0
    return pd.DataFrame({'trial':np.arange(1,len(ts)+1),'input':trainVec[ts-1],'vx':bandVec[trainVec[ts-1]-1]+noiseVec})


def sim_data(dat, c=0.5, lr=0.2, inNodes=7, outNodes=32, trainVec=[5,6,7]):
    global inputNodes, outputNodes
    inputNodes = np.linspace(1, 7, num=inNodes)
    outputNodes = np.linspace(50, 1600, num=outNodes)
    wm = np.zeros((len(outputNodes), len(inputNodes)))
    tt = trainTest_alm(dat, c, lr, wm, trainVec)
    return tt




def RMSE(x,y):
  return np.sqrt(np.mean((x-y)**2))

def RMSE_blocked(x,y,blocks=6):
  df = pd.DataFrame({'x':x,'y':y,'t':np.arange(1,len(x)+1)})
  df['fitBins'] = pd.cut(df['t'],blocks,labels=np.arange(1,blocks+1))
  df = df.groupby('fitBins').agg({'x':np.mean,'y':np.mean}).reset_index()
  return RMSE(df['x'],df['y'])

def MAPE(x, y):
  return np.mean(np.abs((x - y) / y)) * 100

def MedAE(x, y):
  return np.median(np.abs(x - y))

def HuberLoss(x, y, delta = 1):
  error = x - y
  abs_error = np.abs(error)
  loss = np.where(abs_error <= delta, 0.5 * error**2, delta * (abs_error - 0.5 * delta))
  return np.mean(loss)

def sigmoid(x):
  return 1/(1+np.exp(-x))