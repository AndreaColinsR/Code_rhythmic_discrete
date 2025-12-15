function [FRoutput,state]=evalRNN(Input,NetParams)
W=NetParams.W;
O=NetParams.O;
Ob=NetParams.Ob;
B=NetParams.B;
Initial_state=NetParams.S0;

Nunits=size(W,1);
Ntimes=size(Input,1);
NIn=size(Input,2);
FRoutput=zeros(Ntimes,size(O,2));
state=zeros(Nunits,Ntimes);
state(:,1)=Initial_state;
tau=4;


for i=2:Ntimes
    %I=((Input(i,:)+normrnd(0.01,0.01,1,NIn))*B')';
    I=(Input(i,:)*B')';
    state(:,i)=tanh(I+W'*state(:,i-1))./tau;%+normrnd(0,1/1000,Nunits,1);
    %% corrected tanh
    %state(:,i)=(tanh(I+W'*state(:,i-1))+1)./(tau*2);
    FRoutput(i,:)=state(:,i)'*O+Ob;
end
state=state';
end
