import numpy as np
import tensorflow as tf

# Forward simulation of a trained rate-based RNN using an explicit time-step#

def simulateRNN(Bipos,Bdir,B,W,O,Ob,I,tau):
    """
    Parameters
    ----------
    Bipos : array, shape (2, Nunits)
        Input weights for start position (one-hot, 2 dims).
    Bdir : array, shape (2, Nunits)
        Input weights for movement direction (one-hot, 2 dims).
    B : array, shape (1, Nunits)
        Input weights for temporal context / stop-like scalar input.
    W : array, shape (Nunits, Nunits)
        Recurrent weight matrix.
    O : array, shape (Nunits, Noutputs)
        Output weight matrix.
    Ob : array, shape (Noutputs,)
        Output bias.
    I : array, shape (1, Ntimes, Ninputs)
        Input time series for one trial. 
    tau : float
        Time constant.

    Returns
    -------
    Output : array, shape (Ntimes, Noutputs)
        Simulated network output over time.
    state : array, shape (Nunits, Ntimes+1)
        Simulated hidden state trajectory.
    """
    Ntimes = I.shape[1]
    Nunits = W.shape[0]
    Output = np.zeros((Ntimes,O.shape[1]))
    state = np.zeros((Nunits,Ntimes+1))
    #initializer = tf.keras.initializers.Zeros()
    initializer = tf.keras.initializers.RandomNormal(mean=0., stddev=0.1, seed=1)
    S0 = initializer((1, Nunits)).numpy()
    state[:,0:1]=S0.transpose()
    
    
    for t1 in np.arange(0,Ntimes):
        #It=inputs[0,t1,:]
        It=np.matmul(I[0,t1:t1+1,0:2],Bipos[0:2,:])+\
        B[0,:].transpose()*I[0,t1:t1+1,2]+\
        np.matmul(I[0,t1:t1+1,3:5],Bdir[0:2,:])
        
        Rec=np.dot(W.transpose(),state[:,t1])
        ## Output= tanh(I*B+state*W)*O+Ob
        #state[:,t1+1]=(np.tanh(It+Rec)+1)/(tau*2)
        state[:,t1+1]=np.tanh(It+Rec)/tau

        #output = (-prev_output+ops.matmul(prev_output, ops.tanh(self.recurrent_kernel))+ I)/self.tau
        Output[t1,:]=np.dot(state[:,t1+1],O)+Ob

    return Output,state