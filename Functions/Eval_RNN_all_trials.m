function [scores,trials_idx,R2,FRoutput,all_state,Edit_output]=Eval_RNN_all_trials(Input,NetParams,Output,exec,do_plot)
% EVAL_RNN_ALL_TRIALS Evaluates a recurrent neural network (RNN) across all
% trials, computing outputs, hidden states, and performance metrics.
%
% This function simulates the RNN trial by trial, using input signals and
% network parameters, and returns firing rate outputs, hidden states,
% principal component scores, and the trial-aligned outputs. It can
% optionally generate visualisations of inputs, outputs, neural states,
% and principal components.
%
%
%  [scores,trials_idx,R2,FRoutput,all_state,Edit_output] = ...
%      EVAL_RNN_ALL_TRIALS(Input,NetParams,Output,exec,idx_trial_test,do_plot)
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
%   Output          : [T x Ntrials x Nout] array
%       Target outputs corresponding to the input trials.
%
%   exec            : [Ntrials x K] matrix
%       Execution flags or condition descriptors for each trial.
%
%   idx_trial_test  : [Ntrials x 1] vector
%       Condition label or index for each trial, used for rotation of inputs
%       if applicable.
%
%   do_plot         : scalar (0 or 1)
%       Flag indicating whether to generate visualisations of inputs,
%       outputs, neural states, and principal components.
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
%   Edit_output     : [Npoints x Nout] matrix
%       Target outputs reshaped and concatenated across all trials.



Ntimes=size(Input,1);
Ntrials=size(Input,2);
Nunits=size(NetParams.W,1);
all_state=nan(Ntimes*Ntrials,Nunits);
FRoutput=nan(Ntimes*Ntrials,size(Output,3));
Edit_output=nan(Ntimes*Ntrials,size(Output,3));
trials_idx=nan(Ntimes*Ntrials,1);
t=1:Ntimes;

Output=reshape(Output,Ntrials*Ntimes,size(Output,3));

%% Evaluate all trials
counter=1;
%% define initial state as the fixed point 
In=zeros(300,size(Input,3));
[~,all_states_this_trial]=evalRNN(In,NetParams);
NetParams.S0=all_states_this_trial(end,:)';

for itrials=1:Ntrials
    I_this_trial=squeeze(Input(:,itrials,:));
    mov_end=find(exec(itrials,:)>0,1,'last')+40;
   
   [FR_this_trial,all_states_this_trial]=evalRNN(I_this_trial,NetParams);

   FRoutput(counter+1:counter+mov_end,:)=FR_this_trial(1:mov_end,:);
   all_state(counter+1:counter+mov_end,:)=all_states_this_trial(1:mov_end,:);
   trials_idx(counter+1:counter+mov_end)=itrials;

   Edit_output(counter+1:counter+mov_end,:)=Output(Ntimes*(itrials-1)+1:(Ntimes*(itrials-1)+mov_end),:);


   counter=counter+mov_end;

end

% delete all nans
FRoutput(isnan(sum(FRoutput,2)),:)=[];
all_state(isnan(sum(all_state,2)),:)=[];
trials_idx(isnan(sum(trials_idx,2)),:)=[];
Edit_output(isnan(sum(Edit_output,2)),:)=[];


R2=corr(Edit_output(:),FRoutput(:));

[~,scores]=pca(all_state);

if do_plot
    figure
    colour_trial=copper(Ntrials);
    input_names=cell(Ntrials,1);

    for itrials=1:Ntrials
        thistrial=find(trials_idx==itrials);

        subplot(4,2,1)
        hold on
        plot(t,squeeze(Input(:,itrials,:)),'Color',colour_trial(itrials,:))
        title('Inputs')
        input_names{itrials}=['Trial = ' num2str(itrials)];
        

        subplot(4,2,3)
        hold on
        plot(FRoutput(thistrial),'Color',colour_trial(itrials,:))
        title('Outputs')

        if itrials==Ntrials
            subplot(2,2,2)
            plot(all_state(thistrial,:))
            title('Neural states last trial')
            ylabel('FR')
            xlabel('Time')

            subplot(2,2,3)
            hold on
            plot(scores(thistrial,1))
            plot(scores(thistrial,2))
            plot(scores(thistrial,3))

            legend('PC 1','PC 2','PC 3')
            title('PCs last trial')
            xlabel('Time')
        end

        subplot(2,2,4)
        hold on
        plot3(scores(thistrial,1),scores(thistrial,2),scores(thistrial,3),'Color',colour_trial(itrials,:))
        plot3(scores(thistrial(1),1),scores(thistrial(1),2),scores(thistrial(1),3),'o','Color',colour_trial(itrials,:))
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
    end

    subplot(4,2,1)
    legend(input_names)
end



end

