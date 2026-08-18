function Supplementary_variants_models_SMA_same_controls
% SUPPLEMENTARY_SMA_EMG_output
%
% This function evaluates and visualises the temporal context structure
% learned by trained recurrent neural networks (RNNs) for a single animal.
% It compares discrete and rhythmic movement conditions under two training
% regimes:
%       (1) RNN trained to output EMG dynamics with an input indicating
%       the temporal context (same input as SMA same-hyp model)
%       (2) RNN trained to output EMG dynamics with an input indicating
%       the temporal context and movement type
%
%
% No input arguments.
% No output arguments.
%
% Andrea Colins
% 05/06/2026

% select animal and condition to evaluate in both example RNN
animal = {'Cousteau','Drake'};
region_name = 'SMA';

i_dir = 1;
i_pos = 1;


Nfiles = 20;


R2_same = nan(Nfiles,2);
CC_EMG_same = nan(Nfiles,2,2);
CC_Internal_same = nan(Nfiles,2);

R2_different = nan(Nfiles,2);
CC_EMG_different = nan(Nfiles,2,2);
CC_Internal_different = nan(Nfiles,2);


R2 = nan(Nfiles,2);
CC_EMG = nan(Nfiles,2,2);
CC_Internal = nan(Nfiles,2);

R2_constant = nan(Nfiles,2);
CC_EMG_constant = nan(Nfiles,2,2);
CC_Internal_constant = nan(Nfiles,2);


fig1=figure;
do_plot=0;
for i_animal=1:length(animal)
    for i =1:Nfiles

        load(['.\Output_files\scores_' animal{i_animal} '_'  region_name '.mat'],'scores','idx_dir','idx_pos','idx_dist')

        cond_idx_M1=[idx_dir,idx_pos,idx_dist];
        %% Evaluate RNNs from same-control model


        %% Plot selected examples
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %(1) RNN trained to output EMG dynamics (same-control)

        training = 13:20;
        test = 1:12;


        load(['.\Output_files\RNNs_Inputs\SMA_' animal{i_animal} '_same.mat'],'exec','idx_current_cycle')
        
        filename = ['Output_files\TrainedRNNs\SMA_same\Trained_M1_same_' animal{i_animal} '_' num2str(i) '.mat'];

        [R2_same(i,i_animal),CC_EMG_same(i,:,i_animal),CC_Internal_same(i,i_animal)]=plot_dynamics(i,i_dir,i_pos,scores,cond_idx_M1,exec,idx_current_cycle,training,test,filename,do_plot);
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %(2) RNN trained to output EMG dynamics (different-control)
        
        % Define train and test trials 
        training=[0,1,2,3,16,17,18,19]+1;
        test = 5:16;

        load(['.\Output_files\RNNs_Inputs\SMA_' animal{i_animal} '_different.mat'],'exec','idx_current_cycle')
        filename = ['Output_files\TrainedRNNs\SMA_different\Trained_M1_different_' animal{i_animal} '_'  num2str(i) '.mat'];

        [R2_different(i,i_animal),CC_EMG_different(i,:,i_animal),CC_Internal_different(i,i_animal)]=plot_dynamics(i,i_dir,i_pos,scores,cond_idx_M1,exec,idx_current_cycle,training,test,filename,do_plot);


        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (3) RNN trained to output M1 dynamics with an input indicating
        % the temporal context only during movement preparation
     
        % Define train and test trials 
        training = 13:20;
        test = 1:12;

        load(['.\Output_files\RNNs_Inputs\SMA_' animal{i_animal} '_same_peak_onset.mat'],'exec','idx_current_cycle')
        
        filename = ['Output_files\TrainedRNNs\SMA_same_variant_slope\Trained_M1_same_variant_slope_' animal{i_animal} '_' num2str(i) '.mat'];

        [R2(i,i_animal),CC_EMG(i,:,i_animal),CC_Internal(i,i_animal)]=plot_dynamics(i,i_dir,i_pos,scores,cond_idx_M1,exec,idx_current_cycle,training,test,filename,do_plot);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (4) RNN trained to output EMG constant input for rhythmic,
        % different height of discrete 
    
        %[R2(i,i_animal),CC_EMG(i,:,i_animal),CC_Internal(i,i_animal)]=plot_dynamics_different_height(animal{i_animal},i,i_dir,i_pos,scores,cond_idx_M1);
        
        % Define train and test trials 
        training = 13:20;
        test = 1:12;
        load(['.\Output_files\RNNs_Inputs\SMA_' animal{i_animal} '_same_height_slope.mat'],'exec','idx_current_cycle')
        filename = ['Output_files\TrainedRNNs\SMA_same_variant_height\Trained_M1_same_variant_height_' animal{i_animal} '_' num2str(i) '.mat'];

        [R2_constant(i,i_animal),CC_EMG_constant(i,:,i_animal),CC_Internal_constant(i,i_animal)]=plot_dynamics(i,i_dir,i_pos,scores,cond_idx_M1,exec,idx_current_cycle,training,test,filename,do_plot);
         
      
    end
end

threshold=0.8;

successful_same1=R2_same(:,1)>=threshold;
successful_same2=R2_same(:,2)>=threshold;
disp(" Same-control ")
disp(['Animal ' animal{1} ' has ' num2str(nnz(successful_same1)) 'successfully trained networks'])
disp(['Animal ' animal{2} ' has ' num2str(nnz(successful_same2)) 'successfully trained networks'])
disp("")

successful_different1=R2_different(:,1)>=threshold;
successful_different2=R2_different(:,2)>=threshold;
disp(" Different-control")
disp(['Animal ' animal{1} ' has ' num2str(nnz(successful_different1)) 'successfully trained networks'])
disp(['Animal ' animal{2} ' has ' num2str(nnz(successful_different2)) 'successfully trained networks'])
disp("")


successful_1=R2(:,1)>=threshold;
successful_2=R2(:,2)>=threshold;
Nsuc_1=nnz(successful_1);
Nsuc_2=nnz(successful_2);
disp(" Variant same peak onset ")
disp(['Animal ' animal{1} ' has ' num2str(Nsuc_1) 'successfully trained networks'])
disp(['Animal ' animal{2} ' has ' num2str(Nsuc_2) 'successfully trained networks'])
disp("")

successful_constant1=R2_constant(:,1)>=threshold;
successful_constant2=R2_constant(:,2)>=threshold;
disp(" Variant same height slope ")
disp(['Animal ' animal{1} ' has ' num2str(nnz(successful_constant1)) 'successfully trained networks'])
disp(['Animal ' animal{2} ' has ' num2str(nnz(successful_constant2)) 'successfully trained networks'])
disp("")

figure(fig1)
subplot(3,3,1)

plot_fancy_errorbars(0,R2_same(:),[0.9 0.5 0.5])
plot_fancy_errorbars(1,R2_different(:),[0.5 0.9 0.5])
plot_fancy_errorbars(2,R2(:),[0.9 0.5 0.5])
plot_fancy_errorbars(3,R2_constant(:),[0.9 0.5 0.5])

ylim([0,1])
xlim([-0.5 3.5])
xticks([0 1 2 3 4])
xticklabels({'Same','Different','Different height','Same height','Same slope'})
ylabel('Correlation between RNN''s output and M1')
box off

[~,p_different_height] = ttest2(R2_different(:),R2(:));
[~,p_different_peak] = ttest2(R2_different(:),R2_constant(:));

text(2,0.1,num2str(p_different_height))
text(3,0.2,num2str(p_different_peak))

subplot(3,3,2)
title('CCA Output')

CCA_same = [CC_EMG_same(:,2,1);CC_EMG_same(:,2,2)];
CCA_different = [CC_EMG_different(:,2,1);CC_EMG_different(:,2,2)];
CCA_height = [CC_EMG(:,2,1);CC_EMG(:,2,2)];
CCA_peak = [CC_EMG_constant(:,2,1);CC_EMG_constant(:,2,2)];


plot_fancy_errorbars(0,CCA_same,[0.9 0.5 0.5])
plot_fancy_errorbars(1,CCA_different,[0.5 0.9 0.5])
plot_fancy_errorbars(2,CCA_height,[0.9 0.5 0.5])
plot_fancy_errorbars(3,CCA_peak,[0.9 0.5 0.5])

ylim([0,1])
xlim([-0.5 3.5])
xticks([0 1 2 3 4])
xticklabels({'Same','Different','Same height','Same slope'})
box off

[~,p_different_height] = ttest2(CCA_different,CCA_height);
[~,p_different_peak] = ttest2(CCA_different,CCA_peak);

text(2,0.1,num2str(p_different_height))
text(3,0.2,num2str(p_different_peak))

subplot(3,3,3)

title('CCA Internal dynamics')

CCAI_same = [CC_Internal_same(:,1);CC_Internal_same(:,2)];
CCAI_different = [CC_Internal_different(:,1);CC_Internal_different(:,2)];
CCAI_height = [CC_Internal(:,1);CC_Internal(:,2)];
CCA_peak = [CC_Internal_constant(:,1);CC_Internal_constant(:,2)];

plot_fancy_errorbars(0,CCAI_same,[0.9 0.5 0.5])
plot_fancy_errorbars(1,CCAI_different,[0.5 0.9 0.5])
plot_fancy_errorbars(2,CCAI_height ,[0.9 0.5 0.5])
plot_fancy_errorbars(3,CCA_peak,[0.9 0.5 0.5])


ylim([0,1])
xlim([-0.5 4.5])
xticks([0 1 2 3 4])
xticklabels({'Same','Different','Different height','Same slope','Same height'})
ylabel('Correlation between RNN''s output and SMA')
box off

[~,p_different_height] = ttest2(CCAI_different,CCAI_height);
[~,p_different_peak] = ttest2(CCAI_different,CCA_peak);


text(2,0.1,num2str(p_different_height))
text(3,0.2,num2str(p_different_peak))


end


function [R_final,CC_EMG,corr_CC]=plot_dynamics(i,i_dir,i_pos,scores_M1,cond_idx_M1,exec,idx_current_cycle,training,test,filename,do_plot)

%training = 13:20;
%test = 1:12;
exec=[exec(training,:);exec(test,:)];
idx_current_cycle=[idx_current_cycle(training,:);idx_current_cycle(test,:)];

info = load_RNN_info(filename);
[scores,trials_idx,R2,~,Output_edited,Output_RNN] = Eval_RNN_all_conditions(info.Input,info.Output,info.Test_input,info.Test_Outputs,info.net_params,exec,info.idx_conditions_train,info.idx_conditions_test,0);
R_final=R2(2);

[idx_dir,idx_pos,~,idx_dist]=find_idx_conds(trials_idx,info.idx_conds_all,idx_current_cycle);
%% Between EMG and RNN's output
idx_all_cond=[idx_dir,idx_pos,idx_dist];

% train
idx_train=ismember(idx_dist,unique(info.idx_cycles_train));
r_EMG_train=CCA_RNN_M1(Output_RNN(idx_train,:),idx_all_cond(idx_train,:),Output_edited(idx_train,:),idx_all_cond(idx_train,:));
% test
idx_test=ismember(idx_dist,unique(info.idx_cycles_test));
r_EMG_test=CCA_RNN_M1(Output_RNN(idx_test,:),idx_all_cond(idx_test,:),Output_edited(idx_test,:),idx_all_cond(idx_test,:));

CC_EMG=[mean(r_EMG_train),mean(r_EMG_test)];


%% internal dynamics
%% only for test trials
idx_test_M1=ismember(cond_idx_M1(:,3),unique(info.idx_cycles_test));
r=CCA_RNN_M1(scores(idx_test,:),idx_all_cond(idx_test,:),scores_M1(idx_test_M1,:),cond_idx_M1(idx_test_M1,:));
corr_CC=mean(r);

%% plots
if R_final>=0.8 && do_plot
    figure
    plot_dynamics_and_output(idx_dir,i_dir,idx_pos,i_pos,idx_dist,Output_edited,Output_RNN,R_final,i,scores)
    
end

end

function plot_dynamics_and_output(idx_dir,i_dir,idx_pos,i_pos,idx_dist,Output_edited,Output_RNN,R_final,i,scores)
discrete = find(idx_dir == i_dir & idx_pos == i_pos & idx_dist == 0.5);
rhythmic = find(idx_dir == i_dir & idx_pos== i_pos & idx_dist == 7);
Ndist = numel(unique(idx_dist));
colours=plasma(Ndist);


subplot(4,3,7)
hold on
plot(Output_edited(discrete,1),'Color',colours(1,:)*1.8) %% PC1 M1
plot(Output_RNN(discrete,1),'Color',colours(1,:))  %% Output RNN
box off
legend('EMG','RNN output')

subplot(4,3,10)
hold on 
plot(Output_edited(rhythmic,1),'Color',colours(end,:)/2)
plot(Output_RNN(rhythmic,1),'Color',colours(end,:))
box off

legend('M1','RNN output')
title(['R output = ',num2str(R_final)])

subplot(2,3,5)

title(['Recording number = ' num2str(i)])
hold on
plot3(scores(rhythmic(1),1),scores(rhythmic(1),2),scores(rhythmic(1),3),'o','MarkerFaceColor',colours(end,:),'MarkerEdgeColor',colours(end,:))
plot3(scores(rhythmic,1),scores(rhythmic,2),scores(rhythmic,3),'Color',colours(end,:))

plot3(scores(discrete,1),scores(discrete,2),scores(discrete,3),'Color',colours(1,:))
plot3(scores(discrete(1),1),scores(discrete(1),2),scores(discrete(1),3),'o','MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(end,:))

box off

subplot(2,3,6)
plot(scores(discrete,1),'Color',colours(1,:),'LineWidth',3)
hold on
plot(scores(rhythmic,1),'Color',colours(end,:))
title('RNN scores')
box off

end