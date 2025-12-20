function [Output,RNNstates]=evalRNN(Input,NetParams)
% EVALRNN Simulates a recurrent neural network (RNN) given input signals
% and network parameters, returning network outputs and hidden states. This
% function propagates the input signals through the RNN using the
% trained weights and biases, computing the hidden states at each time
% step and the corresponding network outputs. The network uses a
% hyperbolic tangent activation function with a time constant. The time
% constant is fixed at 40 ms, asuming a timebin of 10 ms. 
%
%
%  [Output,RNNstates] = EVALRNN(Input,NetParams)
%
%   INPUTS
%   ------
%   Input       : [T x Nin] matrix
%       Input signals for the RNN, where T is the number of time samples
%       and Nin is the number of input channels.
%
%   NetParams   : structure
%       Structure containing the trained RNN parameters, with fields:
%         - W             : recurrent weight matrix [Nunits x Nunits]
%         - O             : output weight matrix [Nunits x Nout]
%         - Ob            : output bias vector [1 x Nout]
%         - B             : input weight matrix [Nin x Nunits]
%         - Initial_state : initial hidden state vector [Nunits x 1]
%
%   OUTPUTS
%   -------
%   Output      : [T x Nout] matrix
%       Network outputs at each time step, where Nout is the number of
%       output channels.
%
%   RNNstates  : [T x Nunits] matrix
%       Hidden states of the RNN at each time step, with Nunits being the
%       number of recurrent units. Each row corresponds to a time sample.
%
% Andrea Colins
% 19/12/2025



% Get parameters
W=NetParams.W;
O=NetParams.O;
Ob=NetParams.Ob;
B=NetParams.B;
Initial_state=NetParams.S0;

Nunits=size(W,1);
Ntimes=size(Input,1);
Output=zeros(Ntimes,size(O,2));
RNNstates=zeros(Nunits,Ntimes);
RNNstates(:,1)=Initial_state;
tau = 4; % time constant 4 equals 40 ms since dt = 10 ms


for i=2:Ntimes
    
    I=(Input(i,:)*B')';
    
    RNNstates(:,i)=tanh(I+W'*RNNstates(:,i-1))./tau;
   
    Output(i,:)=RNNstates(:,i)'*O+Ob;
end
RNNstates=RNNstates';
end
