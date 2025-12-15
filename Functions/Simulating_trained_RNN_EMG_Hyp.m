function [scores,trials_idx,R2,FRoutput,all_state,Edit_output] = Simulating_trained_RNN_EMG_Hyp(Input,NetParams,Output,exec,idx_trial_test,do_plot)
%% Not a clue why I need this function. Check if I can get rid of it. 

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


Orig_pars=NetParams;

%if size(NetParams.B,2)>6
if size(NetParams.B,2)>10
    rotate=true;
else
    rotate=false;
end

for itrials=1:Ntrials
    I_this_trial=squeeze(Input(:,itrials,:));
    mov_end=find(exec(itrials,:)>0,1,'last')+40;
% 
%     if ~isnan(idx_trial_test(itrials,1)) && idx_trial_test(itrials,1)<2
%          I_this_trial(30:70,7)=0.1;
%      elseif ~isnan(idx_trial_test(itrials,1)) && idx_trial_test(itrials,1)>=2
%         I_this_trial(30:70,6)=0.1;
%      end
    
    if rotate
        NetParams=rotate_angle_for_test(NetParams,idx_trial_test(itrials,1));
    end

   [FR_this_trial,all_states_this_trial]=evalRNN(I_this_trial,NetParams);
   
    
   %% restore original value to rotate the correct angle in the next iteration
   if rotate
    NetParams=Orig_pars;
   end

   FRoutput(counter+1:counter+mov_end,:)=FR_this_trial(1:mov_end,:);
   all_state(counter+1:counter+mov_end,:)=all_states_this_trial(1:mov_end,:);
   trials_idx(counter+1:counter+mov_end)=itrials;

   Edit_output(counter+1:counter+mov_end,:)=Output(Ntimes*(itrials-1)+1:(Ntimes*(itrials-1)+mov_end),:);
   % [FRoutput(Ntimes*(itrials-1)+1:Ntimes*itrials,:),all_state(Ntimes*(itrials-1)+1:Ntimes*itrials,:)]=evalRNN(I_this_trial,NetParams);
   %trials_idx(Ntimes*(itrials-1)+1:Ntimes*itrials)=itrials;

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

function NetParams=rotate_angle_for_test(NetParams,idx_trial_test)

%alpha=pi/2-(pi/2)*([0.5 1 2 4 7]-0.5)/6.5;


alpha=pi/2-(pi/2)./(1+exp(4-2*[0.5 1 2 4 7]));

%alpha=(exp(0.5*(0.5-[0.5 1 2 4 7])))*pi/2;
%Ncycle=[0.5 1 2 4 7];

%% Linear 
%alpha=fliplr(Ncycle);
%% log
%alpha=log(Ncycle(end))-log(Ncycle);
%% decaying exp
%alpha=(exp(1./((Ncycle)))-exp(1./((Ncycle(end)))));


%alpha=(pi/2)*alpha./max(alpha);

if isnan(idx_trial_test) || idx_trial_test==0.5 || idx_trial_test==7
    return
elseif idx_trial_test==1
    NetParams.B(:,7)=rotate_n_dimensional_vector(NetParams.B(:,6),alpha(2));
elseif idx_trial_test==2
    NetParams.B(:,6)=rotate_n_dimensional_vector(NetParams.B(:,6),alpha(3));
elseif idx_trial_test==4
    NetParams.B(:,6)=rotate_n_dimensional_vector(NetParams.B(:,6),alpha(4));

end

end


