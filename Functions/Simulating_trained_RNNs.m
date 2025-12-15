function [scores,trials_idx,R2,FRoutput,all_state,Coutput_reformat] = Simulating_trained_RNNs(Input,NetParams,Cortical_Output,exec)
%% Not a clue why I need this function. Check if I can get rid of it. 

Ntimes=size(Input,1);
Ntrials=size(Input,2);
Nunits=size(NetParams.W,1);
all_state=nan(Ntimes*Ntrials,Nunits);
FRoutput=nan(Ntimes*Ntrials,size(Cortical_Output,3));
Coutput_reformat=nan(Ntimes*Ntrials,size(Cortical_Output,3));
trials_idx=nan(Ntimes*Ntrials,1);

Cortical_Output=reshape(Cortical_Output,Ntrials*Ntimes,size(Cortical_Output,3));

%% Evaluate all trials

%% define initial state as the fixed point 
In=zeros(300,size(Input,3));
[~,RNN_states_this_trial]=evalRNN(In,NetParams);
NetParams.S0=RNN_states_this_trial(end,:)';

counter = 1;
for itrials=1:Ntrials

    % Format input for this condition
    I_this_trial=squeeze(Input(:,itrials,:));
    
    % Evaluate this condition
    [Output_this_trial,RNN_states_this_trial] = evalRNN(I_this_trial,NetParams);
   
   mov_end = find(exec(itrials,:)>0,1,'last') + 40;

   FRoutput(counter+1:counter+mov_end,:)=Output_this_trial(1:mov_end,:);
   all_state(counter+1:counter+mov_end,:)=RNN_states_this_trial(1:mov_end,:);
   trials_idx(counter+1:counter+mov_end)=itrials;
   % change format of cortical activity
   Coutput_reformat(counter+1:counter+mov_end,:) = Cortical_Output(Ntimes*(itrials-1)+1:(Ntimes*(itrials-1)+mov_end),:);

   counter=counter+mov_end;

end

% delete all nans
FRoutput(isnan(sum(FRoutput,2)),:)=[];
all_state(isnan(sum(all_state,2)),:)=[];
trials_idx(isnan(sum(trials_idx,2)),:)=[];
Coutput_reformat(isnan(sum(Coutput_reformat,2)),:)=[];

R2 = corr(Coutput_reformat(:),FRoutput(:));

[~,scores]=pca(all_state);

end