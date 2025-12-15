function [scores_all,trials_idx,R2,states,Output_edit,Output_RNN,r_end]=Eval_RNN_all_conditions(Input,Output,InputTest,OutputTest,NetParams,exec,idx_trial_train,idx_trial_test,do_plot_output)
%% aim: this function evaluates the network using the inputs that were used to trained it and extra inputs to test how well extrapolates to other problems

% idx_trials= Ncycles,Init_pos
idx_all_trials=[idx_trial_train;idx_trial_test];

Ntimes=size(Input,1);
Ntrials=size(idx_all_trials,1);
idc=unique(idx_all_trials(:,1),'rows','stable');
idc_sorted=unique(idx_all_trials(:,1));
Ncycles=size(idc,1);
Noutputs=size(OutputTest,3);



do_plot=0;


%% Evaluate training
[~,trials_idx,R2,Output_RNN,states_train,Output_edit]=Simulating_trained_RNN_EMG_Hyp(Input,NetParams,Output,exec(1:size(idx_trial_train,1),:),idx_trial_train,do_plot);

%% Evaluate testing


[~,trials_idx_test,R2_test,Output_RNNTest,states_test,Output_edittest]=Simulating_trained_RNN_EMG_Hyp(InputTest,NetParams,OutputTest,exec(size(idx_trial_train,1)+1:end,:),idx_trial_test,do_plot);
trials_idx=[trials_idx;trials_idx_test+max(trials_idx)];
Input=[Input,InputTest];
Output=[Output,OutputTest];
Output_edit=[Output_edit;Output_edittest];
Output_RNN=[Output_RNN;Output_RNNTest];
states=[states_train;states_test];
[~,scores_all]=pca(states);


r_end=test_end_time(Input,scores_all,trials_idx);

R2=[R2 R2_test];
%%%%% plot?
%do_plot_output=1;

if do_plot_output
t=(1:Ntimes)*10-1000;
figure
colour_trial=plasma(Ncycles);
colour_init_pos=jet(5);
colour_dir=[1 0 0;0 0 0];



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
        title(['R = ' num2str(R2_trial)])
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


        subplot(3,3,3)
        hold on
        plot3(scores_all(thistrial,1),scores_all(thistrial,2),scores_all(thistrial,3),'Color',colour_trial(i_cycle_sorted,:),'LineWidth',i_initpos^2)
        plot3(scores_all(thistrial(1),1),scores_all(thistrial(1),2),scores_all(thistrial(1),3),'o','Color',colour_trial(i_cycle_sorted,:))
        title(['N cycle = ' num2str(idx_all_trials(itrial,1))])

        subplot(3,3,6)
        hold on
        plot3(scores_all(thistrial,1),scores_all(thistrial,2),scores_all(thistrial,3),'Color',colour_init_pos(idx_all_trials(itrial,2),:),'LineWidth',i_initpos^2)
        plot3(scores_all(thistrial(1),1),scores_all(thistrial(1),2),scores_all(thistrial(1),3),'o','Color',colour_init_pos(idx_all_trials(itrial,2),:))
        title(['N cycle = ' num2str(idx_all_trials(itrial,1))])
        if idx_all_trials(itrial,1)==0.5
        
        subplot(3,3,9)
        hold on
        %plot3(scores_all(thistrial,1),scores_all(thistrial,2),scores_all(thistrial,3),'Color',colour_dir(idx_all_trials(itrial,3),:),'LineWidth',i_initpos^2)
        %plot3(scores_all(thistrial(150:300),1),scores_all(thistrial(150:300),2),scores_all(thistrial(150:300),3),'Color',colour_dir(idx_all_trials(itrial,3),:),'LineWidth',i_initpos^2)
       
        plot3(scores_all(thistrial,1),scores_all(thistrial,2),scores_all(thistrial,3),'Color',colour_dir(idx_all_trials(itrial,3),:),'LineWidth',i_initpos^2)
        %plot3(scores_all(thistrial(1),1),scores_all(thistrial(1),2),scores_all(thistrial(1),3),'o','Color',colour_dir(idx_all_trials(itrial,3),:))
        title(['N cycle = ' num2str(idx_all_trials(itrial,1))])
        end
        end
end
end

end

function mean_r=test_end_time(Input,scores_all,trials_idx)
r=nan(20,1);

N_in=3;


for i=1:20
  
if sum(Input(:,i,N_in))==0
   if any(diff(Input(:,i,6))<0 & diff(Input(:,i,6))>-1)
       N_in=6;
   else
       N_in=7;
   end
end

t_end=find(Input(:,i,N_in)>0,1,'last');
X=scores_all(trials_idx==i,1:10);

[Nt,~]=size(X);
I_end=zeros(Nt,1);

I_end(t_end-30:t_end)=1;

beta = mvregress(X,I_end);

Y=X*beta;
r(i)=corr(Y,I_end);
end

mean_r=mean(r);

end