function [Output,RNNstates]=evalRNN(Input,NetParams)
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
    %I=((Input(i,:)+normrnd(0.01,0.01,1,NIn))*B')';
    I=(Input(i,:)*B')';
    RNNstates(:,i)=tanh(I+W'*RNNstates(:,i-1))./tau;%+normrnd(0,1/1000,Nunits,1);
    %% corrected tanh
    %state(:,i)=(tanh(I+W'*state(:,i-1))+1)./(tau*2);
    Output(i,:)=RNNstates(:,i)'*O+Ob;
end
RNNstates=RNNstates';
end
