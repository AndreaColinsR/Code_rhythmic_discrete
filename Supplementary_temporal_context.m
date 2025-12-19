function Supplementary_temporal_context
% SUPPLEMENTARY_TEMPORAL_CONTEXT
%
% This function evaluates and visualises the temporal context structure
% learned by trained recurrent neural networks (RNNs) for a single animal.
% It compares discrete and rhythmic movement conditions under two training
% regimes: 
%       (1) RNN trained to output M1 dynamics without an input indicating
%       the temporal context 
%       (2) RNN trained to output M1 dynamics with an input indicating
%       the temporal context only during movement preparation
%
%
% No input arguments.
% No output arguments.
%
% Andrea Colins
% 19/12/2025

% select animal and condition to evaluate in both example RNN
animal = 'Drake';
i_dir = 1;
i_pos = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1) RNN trained to output M1 dynamics without an input indicating
%the temporal context
load(['.\Output_files\RNNs_Inputs\M1_' animal '_same.mat'],'exec','idx_current_cycle')

training = 13:20;
test = 1:12;
exec=[exec(training,:);exec(test,:)];
idx_current_cycle=[idx_current_cycle(training,:);idx_current_cycle(test,:)];

filename = 'Output_files\TrainedRNNs\SMA_Control_same\Trained_M1_Hyp_continuousDrake_7.mat';
info = load_RNN_info(filename);
% Evaluate 
[scores,trials_idx] = Eval_RNN_all_conditions(info.Input,info.Output,info.Test_input,info.Test_Outputs,info.net_params,exec,info.idx_conditions_train,info.idx_conditions_test,0);

[idx_dir,idx_pos,~,idx_dist]=find_idx_conds(trials_idx,info.idx_conds_all,idx_current_cycle);

% plot
Ndist = numel(unique(idx_dist));
colours=plasma(Ndist);

subplot(2,3,3)
discrete = find(idx_dir == i_dir & idx_pos == i_pos & idx_dist == 0.5);
rhythmic = find(idx_dir == i_dir & idx_pos== i_pos & idx_dist == 7);

hold on
plot3(scores(rhythmic(1),1),scores(rhythmic(1),2),scores(rhythmic(1),3),'o','MarkerFaceColor',colours(end,:),'MarkerEdgeColor',colours(end,:))
plot3(scores(rhythmic,1),scores(rhythmic,2),scores(rhythmic,3),'Color',colours(end,:))

plot3(scores(discrete,1),scores(discrete,2),scores(discrete,3),'Color',colours(1,:))
plot3(scores(discrete(1),1),scores(discrete(1),2),scores(discrete(1),3),'o','MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(end,:))


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) RNN trained to output M1 dynamics with an input indicating
% the temporal context only during movement preparation

load(['.\Output_files\RNNs_Inputs\M1_' animal '_different.mat'],'exec','idx_current_cycle')
training=[0,1,2,3,16,17,18,19]+1;
test = 5:16 ;
exec=[exec(training,:);exec(test,:)];
idx_current_cycle=[idx_current_cycle(training,:);idx_current_cycle(test,:)];

filename = 'Output_files\TrainedRNNs\SMA_Control_different\Trained_M1_Hyp_separateDrake_4.mat';
info = load_RNN_info(filename);
[scores,trials_idx] = Eval_RNN_all_conditions(info.Input,info.Output,info.Test_input,info.Test_Outputs,info.net_params,exec,info.idx_conditions_train,info.idx_conditions_test,0);

[idx_dir,idx_pos,~,idx_dist]=find_idx_conds(trials_idx,info.idx_conds_all,idx_current_cycle);


subplot(2,3,6)
discrete = find(idx_dir == i_dir & idx_pos == i_pos & idx_dist == 0.5);
rhythmic = find(idx_dir == i_dir & idx_pos== i_pos & idx_dist == 7);

hold on
plot3(scores(rhythmic(1),1),scores(rhythmic(1),2),scores(rhythmic(1),3),'o','MarkerFaceColor',colours(end,:),'MarkerEdgeColor',colours(end,:))
plot3(scores(rhythmic,1),scores(rhythmic,2),scores(rhythmic,3),'Color',colours(end,:))

plot3(scores(discrete,1),scores(discrete,2),scores(discrete,3),'Color',colours(1,:))
plot3(scores(discrete(1),1),scores(discrete(1),2),scores(discrete(1),3),'o','MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(end,:))


end