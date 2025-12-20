function [scores_all,trials_idx,R2,states,Output_edit,Output_RNN]=Eval_RNN_all_conditions(Input,Output,InputTest,OutputTest,NetParams,exec,idx_trial_train,idx_trial_test,do_plot_output)
% EVAL_RNN_ALL_CONDITIONS Evaluates a trained recurrent neural network (RNN)
% across training and testing trials, quantifying performance,
% population dynamics, and generalisation to testing trials.
%
% The function simulates the RNN using both the training inputs and
% additional test inputs, computes performance metrics, extracts network
% states, and performs dimensionality reduction on the population activity.
% Optional visualisations compare target outputs and RNN predictions across
% task conditions.
%
%
%  EVAL_RNN_ALL_CONDITIONS(Input,Output,InputTest,OutputTest,NetParams, ...
%                          exec,idx_trial_train,idx_trial_test,do_plot_output)
%
%   INPUTS
%   ------
%   Input           : [T x Ntrain x Nin] array
%       Input signals used during network training, where T is the number of
%       time samples, Ntrain is the number of training trials, and Nin is the
%       number of input channels.
%
%   Output          : [T x Ntrain x Nout] array
%       Target output signals corresponding to the training inputs.
%
%   InputTest       : [T x Ntest x Nin] array
%       Input signals used to evaluate network generalisation on test trials.
%
%   OutputTest      : [T x Ntest x Nout] array
%       Target output signals corresponding to the test inputs.
%
%   NetParams       : structure
%       Structure containing the trained RNN parameters, with fields:
%         - W             : recurrent weight matrix
%         - O             : output weight matrix
%         - Ob            : output bias vector
%         - B             : input bias vector
%         - Initial_state : initial hidden state of the RNN
%
%   exec            : [Ntrain+Ntest x T] matrix
%       Execution of movement flag. exec equals 1 if animals is 
%       executing a movements, 0 otherwise.
%
%   idx_trial_train : [Ntrain x 3] matrix
%       Condition labels for training trials (e.g. number of cycles,
%       initial position, and movement direction).
%
%   idx_trial_test  : [Ntest x 3] matrix
%       Condition labels for test trials.
%
%   do_plot_output  : scalar (0 or 1)
%       Flag indicating whether to generate output and input visualisations.
%
%   OUTPUTS
%   -------
%   scores_all      : matrix
%       Principal component scores obtained from PCA applied to the RNN
%       population states across all trials.
%
%   trials_idx      : vector
%       Index mapping each time point to its corresponding trial.
%
%   R2              : vector
%       Performance metric (correlation coefficient) for each trial,
%       including both training and testing conditions.
%
%   states          : matrix
%       Concatenated RNN hidden states across all trials and time points.
%
%   Output_edit     : matrix
%       Target outputs reshaped and concatenated across trials.
%
%   Output_RNN      : matrix
%       RNN-generated outputs concatenated across training and testing trials.
%
%
% Andrea Colins Rodriguez
% 19/12/2025

idx_all_trials=[idx_trial_train;idx_trial_test];

Ntimes=size(Input,1);
Ntrials=size(idx_all_trials,1);
idc=unique(idx_all_trials(:,1),'rows','stable');
idc_sorted=unique(idx_all_trials(:,1));
Ncycles=size(idc,1);
Noutputs=size(OutputTest,3);


%% Evaluate training
[~,trials_idx,R2,Output_RNN,states_train,Output_edit] = Simulating_trained_RNNs(Input,NetParams,Output,exec(1:size(idx_trial_train,1),:));

%% Evaluate testing
[~,trials_idx_test,R2_test,Output_RNNTest,states_test,Output_edittest] = Simulating_trained_RNNs(InputTest,NetParams,OutputTest,exec(size(idx_trial_train,1)+1:end,:));


trials_idx=[trials_idx;trials_idx_test+max(trials_idx)];
Input=[Input,InputTest];
Output_edit=[Output_edit;Output_edittest];
Output_RNN=[Output_RNN;Output_RNNTest];
states=[states_train;states_test];


[~,scores_all]=pca(states);

R2=[R2 R2_test];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_plot_output

t=(1:Ntimes)*10-1000;
colour_trial=plasma(Ncycles);

for itrial=1:Ntrials
    
        thistrial=find(trials_idx==itrial);
        npoints=numel(thistrial);
        i_cycle=find(idx_all_trials(itrial,1)==idc);
        i_cycle_sorted=find(idx_all_trials(itrial,1)==idc_sorted); %for colour code
        
        i_initpos=idx_all_trials(itrial,2);
        i_dir=round(idx_all_trials(itrial,3));

        Output_trial=Output_edit(thistrial,:);
        Input_trial=squeeze(Input(:,itrial,:));        
        
        if i_dir==1 && i_initpos==1
        %% plot input
        subplot(Ncycles,3,3*i_cycle-2)
        hold on
        plot(t(1:npoints),Input_trial(1:npoints,1)+4,'Color',colour_trial(i_cycle_sorted,:))
        plot(t(1:npoints),Input_trial(1:npoints,2)+2,'Color',colour_trial(i_cycle_sorted,:),'LineWidth',2)
        plot(t(1:npoints),Input_trial(1:npoints,3),'Color',colour_trial(i_cycle_sorted,:),'LineWidth',2)
        xlim([-1000 4000])
        xlabel('Time from mov onset [ms]')
        R2_trial=corr(reshape(Output_RNN(thistrial,:),npoints*Noutputs,1),Output_trial(:));
        title(['R = ' num2str(round(R2_trial,1))])
        %% plot output
        
        subplot(Ncycles,3,3*i_cycle-1)
        hold on
        
        %% target
        plot(t(1:npoints),Output_trial(:,1)+1.5,'LineWidth',3,'Color',[0.5 0.5 0.5])
        plot(t(1:npoints),Output_trial(:,2)+1,'LineWidth',3,'Color',[0.5 0.5 0.5])
        plot(t(1:npoints),Output_trial(:,3)+0.5,'LineWidth',3,'Color',[0.5 0.5 0.5])
        plot(t(1:npoints),Output_trial(:,4),'LineWidth',3,'Color',[0.5 0.5 0.5])

        %% RNN solution
        plot(t(1:npoints),Output_RNN(thistrial,1)+1.5,'Color',colour_trial(i_cycle_sorted,:))
        plot(t(1:npoints),Output_RNN(thistrial,2)+1,'Color',colour_trial(i_cycle_sorted,:))
        plot(t(1:npoints),Output_RNN(thistrial,3)+0.5,'Color',colour_trial(i_cycle_sorted,:))
        plot(t(1:npoints),Output_RNN(thistrial,4),'Color',colour_trial(i_cycle_sorted,:))
        
        xlim([-1000 4000])
        xlabel('Time from mov onset [ms]')

        end
end
end

end