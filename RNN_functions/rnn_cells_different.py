import tensorflow as tf
from tensorflow import keras
    
# Let's define a RNN Cell, as a layer subclass.

class MinimalRNNCell(keras.layers.Layer):

    def __init__(self, units, tau,Btask,**kwargs):
        
        """
        Parameters
        ----------
        units : int
            Number of recurrent units
        tau : float
            Time constant
        Btask : shape (units, n_tasks)
            Fixed projection matrix for the task (movement-type) input.
        """     

        super().__init__(**kwargs)
        self.units = units
        self.state_size = units
        self.tau = tau
        #self.Bipos = Bipos
        #self.Bdir = Bdir
        self.Btask = Btask

    def build(self, input_shape):
        """
        Create and initialise trainable weights of the RNN cell

        Parameters
        ----------
        input_shape :
        Shape of the input tensor (batch, time, Ninputs)
        In this notebook, Ninputs = 5:
            - 2 inputs encoding start position
            - 1 temporal context / stop input
            - 2 inputs encoding movement direction
        """
         
        Initializer=tf.keras.initializers.RandomNormal(mean=0.0, stddev=0.3, seed=None)
        self.B = self.add_weight(shape=(1, self.units),
                                      initializer=Initializer,
                                      name='Input_weights')
        self.Bdir = self.add_weight(shape=(2, self.units),
                                      initializer='RandomUniform',
                                      name='Input_direction')
        self.Bipos = self.add_weight(shape=(2, self.units),
                                      initializer='RandomUniform',
                                      name='Input_pos')
        self.recurrent_kernel = self.add_weight(
            shape=(self.units, self.units),
            initializer=Initializer,
            name='recurrent_kernel')
        self.built = True

    def call(self, inputs, states):
        """
        Perform one recurrent update of the RNN cell

        Parameters
        ----------
        inputs : tensor, shape (batch, Ninputs)
        Inputs at the current time step
        Expected Ninputs = 7 (start position, temporal context, direction)
            - inputs[:, 0:2] : start position (2D)
            - inputs[:, 2:3] : temporal context / stop signal (1D)
            - inputs[:, 3:5] : movement direction (2D)
            - inputs[:, 5:7] : task / movement-type input (2D, e.g. discrete vs rhythmic)
        states : list of tensor
        Previous hidden state of the network, shape (batch, units)

        Returns
        -------
        output : tensor, shape (batch, units)
        Updated hidden state after one time step
        new_states : list
        List containing the updated hidden state
        """
        prev_output = states[0]
        I = tf.matmul(inputs[:,0:2], self.Bipos[0:2,:])+\
            tf.matmul(inputs[:,2:3], self.B)+\
            tf.matmul(inputs[:,3:5], self.Bdir[0:2,:])+\
            tf.matmul(inputs[:,5:7], self.Btask[0:2,:])

        ## From Russo 2020 
        output = tf.tanh(I + tf.matmul(prev_output, self.recurrent_kernel))/self.tau
        #output = (tf.tanh(I + tf.matmul(prev_output, self.recurrent_kernel))+1)/(self.tau*2)
        return output, [output]

    def get_initial_state(self, inputs=None, batch_size=None, dtype=None):
        """
        Initialise the hidden state of the RNN

        The initial state is drawn from a fixed random distribution to ensure
        reproducibility across training runs.
        """
              
        #initializer = tf.keras.initializers.Zeros()
        initializer = tf.keras.initializers.RandomNormal(mean=0., stddev=0.1, seed=1)
        S0= tf.repeat(initializer((1, self.units)),batch_size,axis=0)
        return S0