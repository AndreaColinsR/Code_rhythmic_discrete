function [scores,trials_idx,R2,FRoutput,all_state,Coutput_reformat] = Simulating_trained_RNNs(Input,NetParams,Cortical_Output,exec)
% SIMULATING_TRAINED_RNNS Evaluates a recurrent neural network (RNN) across all
% trials, computing outputs, hidden states, and performance metrics.
%
% This function simulates the RNN trial by trial, using input signals and
% network parameters, and returns firing rate outputs, hidden states,
% principal component scores, and the trial-aligned outputs. It can
% optionally generate visualisations of inputs, outputs, neural states,
% and principal components.
%
%
%  [scores,trials_idx,R2,FRoutput,all_state,Coutput_reformat] = ...
%      EVAL_RNN_ALL_TRIALS(Input,NetParams,Cortical_Output,exec)
%
%   INPUTS
%   ------
%   Input           : [T x Ntrials x Nin] array
%       Input signals for all trials, where T is the number of time samples,
%       Ntrials is the number of trials, and Nin is the number of input channels.
%
%   NetParams       : structure
%       Structure containing the trained RNN parameters, with fields:
%         - W             : recurrent weight matrix [Nunits x Nunits]
%         - O             : output weight matrix [Nunits x Nout]
%         - Ob            : output bias vector [1 x Nout]
%         - B             : input weight matrix [Nin x Nunits]
%         - Initial_state : initial hidden state vector [Nunits x 1]
%
%   Cortical_Output : [T x Ntrials x Nout] array
%       Target outputs corresponding to the input trials.
%
%   exec            : [Ntrials x T] matrix
%       Execution of movement flag. exec equals 1 if animals is 
%       executing a movements, 0 otherwise.
%
%   OUTPUTS
%   -------
%   scores          : matrix
%       Principal component scores obtained from PCA applied to the
%       concatenated RNN states across all trials.
%
%   trials_idx      : vector
%       Index mapping each time point to its corresponding trial.
%
%   R2              : scalar
%       Overall correlation coefficient between target outputs and RNN outputs.
%
%   FRoutput        : [Npoints x Nout] matrix
%       RNN-generated firing rate outputs concatenated across all trials
%       (Npoints = total time points across all trials).
%
%   all_state       : [Npoints x Nunits] matrix
%       Hidden states of the RNN concatenated across all trials.
%
%   Coutput_reformat     : [Npoints x Nout] matrix
%       Target outputs reshaped and concatenated across all trials.
%
% Andrea Colins Rodriguez
% 19/12/2025


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