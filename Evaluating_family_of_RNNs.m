function [corr_CC,Angle_rotRNN,OutW,MinDist,Dist_prep_onset]=Evaluating_family_of_RNNs(hyp,region_name,my_dir,plot_supp_figs,varargin)

if nargin==6
    ifamily=varargin{2};
    figW=varargin{1};
elseif nargin==5
    ifamily=1;
    figW=varargin{1};
else
    figW=figure;
    ifamily=1;
end

plot_supp = plot_supp_figs.do_plot;

hold on

Nnets=20;
animal = {'Cousteau','Drake'}; % 'E' or 'F'
Nanimal=numel(animal);

SuccesfulNets=zeros(Nanimal,1);
corr_CC=nan(Nnets*Nanimal,7);
R2=nan(Nnets*Nanimal,2);
Mean_OutW = nan(Nnets,2);
Angle_rot = nan(Nnets,2);
Angle_disc_rhyth = nan(Nnets,2);
MinDist=nan(Nnets*Nanimal,2);
Init_cond_t = nan(Nnets*Nanimal*4,5);
Nel=600;
Dist2Att_all=nan(Nel,Nnets*Nanimal,5);
Dist_all_prep=nan(100,3,Nnets*Nanimal);


%% figures format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colour_dist=plasma(5);
if ifamily==1
    colour_animal=[[1 0 0];[1 0.5 0.5]];

elseif ifamily==2
    colour_animal=[[0 1 0];[0.5 1 0.5]];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Evaluate each animal recording
for i_animal=1:Nanimal
    load([ my_dir '\Output_files\scores_' animal{i_animal} '_'  region_name '.mat'],'scores','idx_dir','idx_pos','idx_dist')

    cond_idx_M1=[idx_dir,idx_pos,idx_dist];

    Results = Evaluate_this_animal(animal{i_animal},hyp,figW,ifamily,scores,cond_idx_M1,region_name,plot_supp_figs);

    Angle_rot(:,i_animal)=Results.Angle_rot;
    Angle_disc_rhyth(:,i_animal)=Results.Angle_disc_rhyth;
    SuccesfulNets(i_animal)=Results.SuccesfulNets;
    corr_CC((i_animal-1)*Nnets+1:i_animal*Nnets,:)=Results.corr_CC;
    Percentage_all_prep((i_animal-1)*Nnets+1:i_animal*Nnets,:)=Results.Percentage_all_prep;
    Mean_OutW(:,i_animal)=Results.Mean_OutW;
    Init_cond_t((i_animal-1)*Nnets*4+1:i_animal*Nnets*4,:)=Results.Init_cond_t;
    Dist2Att_all(1:size(Results.Dist2Att_all,1),(i_animal-1)*Nnets+1:i_animal*Nnets,:)=Results.Dist2Att_all;
    Dist_all_prep(:,:,(i_animal-1)*Nnets+1:i_animal*Nnets)=Results.Dist_all_prep;
    R2((i_animal-1)*Nnets+1:i_animal*Nnets,:)=Results.R2;

end

Angle_rotRNN=Angle_rot(:);
OutW=Mean_OutW(:);
figure(figW)

%% CCA Output RNN vs Output cortical region (e.g EMG vs output of RNN corresponding to M1)
subplot(4,4,4)
hold on
plot_fancy_errorbars(ifamily-0.2,corr_CC(:,end-1),colour_animal(1,:))
plot_fancy_errorbars(ifamily+0.2,corr_CC(:,end),colour_animal(1,:))


%% successful nets
% subplot(5,4,8)
% hold on
% plot(ifamily,sum(SuccesfulNets),'.','Color',colour_animal(1,:))

%% CCA RNN Internal dynamics vs cortical region dynamics (e.g M1 vs RNN corresponding to M1)
subplot(4,4,8)
hold on
plot_fancy_errorbars(ifamily,corr_CC(:,1),colour_animal(1,:)) % main number
plot_fancy_errorbars(3,corr_CC(:,end-2),[0 0 0]) %shuffle test
[~,pval_shuffle]=ttest2(corr_CC(:,1),corr_CC(:,end-2));
text(0.8,1-ifamily*0.1,num2str(pval_shuffle,'%.2e'),'Units','normalized')
plot_fancy_errorbars(ifamily+0.5,R2(:,2),colour_animal(1,:))

%% var explained
%subplot(4,4,12)
%hold on
%plot_fancy_errorbars(ifamily,Percentage_all_prep(:,1),colour_animal(1,:))


corr_CC=corr_CC(:,1:end-3);

%% angle between rhythmic activity of different directions
% subplot(5,4,12)
% hold on
% plot_fancy_errorbars(ifamily,Angle_rot(:),colour_animal(1,:))

Dist_prep_onset=squeeze(Dist_all_prep(end,1,:));

%% distaces between trajectories (prep)
subplot(4,4,12+ifamily)
t=fliplr((1:size(Dist_all_prep,1))*-10);
hold on
plot_std(t,mean(Dist_all_prep(:,1,:),3,'omitnan'),std(Dist_all_prep(:,1,:),[],3,'omitnan'),[0.1 0.1 0.1])%% N cycle
plot_std(t,mean(Dist_all_prep(:,2,:),3,'omitnan'),std(Dist_all_prep(:,2,:),[],3,'omitnan'),[0.5 0.6 0.9])%% pos
plot_std(t,mean(Dist_all_prep(:,3,:),3,'omitnan'),std(Dist_all_prep(:,3,:),[],3,'omitnan'),[0.5 0.5 0.5])%% dir
ylim([-0.1 1.2])

if plot_supp


    %% Angle between rotations in discrete and rhythmic
    if strcmp(region_name,'M1')
    figure(plot_supp_figs.LDS)
    subplot(2,4,8)
    hold on
    plot_fancy_errorbars(ifamily,Angle_disc_rhyth(:),colour_animal(1,:))

    end
    


    % var explained all pars
    if strcmp(region_name,'M1')
        figure(plot_supp_figs.dPCA_M1)
    else
        figure(plot_supp_figs.dPCA_SMA)
    end

    subplot(5,3,12+ifamily)
    hold on
    errorbar(1:4,mean(Percentage_all_prep,'omitnan'),std(Percentage_all_prep,[],'omitnan'),'.-','Color',colour_animal(1,:))
    xticks(1:4)
    xticklabels({'N cycle','Direction','Pos','CI'})
    figure(figW)
end


if strcmp(region_name,'SMA')

    %% SMA plots
    subplot(4,4,12)
    hold on
    errorbar(mean(Init_cond_t,'omitnan'),[0.5 1 2 4 7],std(Init_cond_t,'omitnan'),'.-','Color',colour_animal(1,:))
    xlim([0 1])
    yticks([0.5 1 2 4 7])
    ylabel('Number of cycles')
    xlabel('Position trajectory at mov onset')
    xlim([-0.1 1])

else
    %% M1 plots

    % distace to the attractor
    subplot(4,4,8+ifamily)
    hold on
    t2=(0:Nel-1)*10-1000;
    for i_dist=1:5
        plot(t2,mean(Dist2Att_all(:,:,i_dist),2,'omitnan'),'Color',colour_dist(i_dist,:))
    end
    tmpminD=min(Dist2Att_all(:,:,1:2));
    MinDist= squeeze(tmpminD);
end
end



function Results = Evaluate_this_animal(animal,hyp,figW,ifamily,scores_M1,cond_idx_M1,region_name,plot_supp_figs)

dists=[0.5 1 2 4 7];
Nel=600;
start=1;
ndims=4;
Npoints=[1 10];

ff=dir(['.\*' hyp '' animal '*']);
Nnetworks=size(ff,1);

figure(figW)

if strcmp(hyp,'separate')
    load(['.\..\..\RNNs_Inputs\M1_' animal '_different.mat'],'exec','idx_current_cycle')
    training=[0,1,2,3,16,17,18,19]+1;
    test=5:16;
    plot_column  = 2;

else
    load(['.\..\..\RNNs_Inputs\M1_' animal '_same.mat'],'exec','idx_current_cycle')
    training=13:20;
    test=1:12;

    plot_column = 1;

    % for inverse case
    %training=1:12;
    %test=13:20;
end


exec=[exec(training,:);exec(test,:)];
idx_current_cycle=[idx_current_cycle(training,:);idx_current_cycle(test,:)];
Dist_all_prep=nan(100,3,Nnetworks);
Percentage_all_prep=nan(Nnetworks,4);
R2=nan(Nnetworks,2);

corr_CC=nan(Nnetworks,1);
corr_CC_control=nan(Nnetworks,1);
CC_EMG=nan(Nnetworks,2);
rprep=nan(Nnetworks,1);
rexec=nan(Nnetworks,1);
rend=nan(Nnetworks,1);
Angle_rot=nan(Nnetworks,1);
Angle_disc_rhyth=nan(Nnetworks,1);
Init_cond_t=nan(Nnetworks*4,5);
Dist2Att_all=nan(Nel,20,numel(dists));

plot_supp = plot_supp_figs.do_plot;

if plot_supp_figs.do_plot == 1
    if strcmp(region_name,'M1')
        prep_fig = plot_supp_figs.Prep_M1;
        dpca_fig = plot_supp_figs.dPCA_M1;
    else
        prep_fig = plot_supp_figs.Prep_SMA;
        dpca_fig = plot_supp_figs.dPCA_SMA;
    end
else
    prep_fig = [];
    dpca_fig = [];
end


for iNet=1:Nnetworks
    if iNet==1  && plot_supp && strcmp(animal,'Cousteau') && strcmp(region_name,'M1')
        figure(plot_supp_figs.Ouputs{ifamily})
        do_plot_output=1;
    else
        do_plot_output=0;
    end

    ff(iNet).name

    load(ff(iNet).name,'inputs','B','W','O','Ob','Output','Test_input','Test_Outputs','idx_cycles_train','idx_pos_train','idx_cycles_test','idx_pos_test','Bipos','Bdir','Btask','Bend','Initial_state','idx_dir_test','idx_dir_train')
    Input = permute(inputs,[2 1 3]);
    Output = permute(Output,[2 1 3]);
    Test_input = permute(Test_input,[2 1 3]);
    Test_Outputs = permute(Test_Outputs,[2 1 3]);

    idx_pos_test=double(idx_pos_test);
    idx_dir_test=double(idx_dir_test);
    idx_pos_train=double(idx_pos_train);
    idx_dir_train=double(idx_dir_train);

    if exist('Btask','var')
        if ~exist('Bend','var')
            %Bend=[B;B];
            net_params.B=[Bipos',B',Bdir',Btask'];
        else
            net_params.B=[Bipos',Bend',Bdir',Btask'];
        end
    else
        net_params.B=[Bipos',B',Bdir'];
    end
    net_params.S0=Initial_state';
    net_params.W=W;
    net_params.O=O;
    net_params.Ob=Ob;

    Mean_OutW(iNet)=norm(O);
    idx_conditions_train=[idx_cycles_train,idx_pos_train,idx_dir_train];
    idx_conditions_test=[idx_cycles_test,idx_pos_test,idx_dir_test];
    idx_conds_all=[idx_conditions_train;idx_conditions_test];

    [scores,trials_idx,R2(iNet,:),states,Output_edited,Output_RNN] = Eval_RNN_all_conditions(Input,Output,Test_input,Test_Outputs,net_params,exec,idx_conditions_train,idx_conditions_test,do_plot_output);


    % post process only if succesful
    if R2(iNet,2)>=0.8

        [idx_dir,idx_pos,idx_Ncycle,idx_dist,~,exec2]=find_idx_conds(trials_idx,idx_conds_all,idx_current_cycle);

        %% LDS
        if strcmp(region_name,'M1') && plot_supp
            [epochs_all(:,:,iNet),t_end]=LDS_RNN(states,idx_dir,idx_pos,idx_dist);
        end

        %% CCA

        %% Between EMGs train
        idx_all_cond=[idx_dir,idx_pos,idx_dist];
        idx_train=ismember(idx_dist,unique(idx_cycles_train));
        r_EMG_train=CCA_RNN_M1(Output_RNN(idx_train,:),idx_all_cond(idx_train,:),Output_edited(idx_train,:),idx_all_cond(idx_train,:));

        %% Between EMG test
        idx_test=ismember(idx_dist,unique(idx_cycles_test));
        r_EMG_test=CCA_RNN_M1(Output_RNN(idx_test,:),idx_all_cond(idx_test,:),Output_edited(idx_test,:),idx_all_cond(idx_test,:));

        CC_EMG(iNet,:)=[mean(r_EMG_train),mean(r_EMG_test)];


        %% control to show that a random model (shuffled RNN PCs) is not
        % well correlated
        r_control=CCA_RNN_M1(scores(:,randperm(size(W,1))),idx_all_cond,scores_M1,cond_idx_M1);
        corr_CC_control(iNet)=mean(r_control);

        %r=CCA_RNN_M1(scores,idx_all_cond,scores_M1,cond_idx_M1);
        %% only for test trials
        idx_test_M1=ismember(cond_idx_M1(:,3),unique(idx_cycles_test));
        r=CCA_RNN_M1(scores(idx_test,:),idx_all_cond(idx_test,:),scores_M1(idx_test_M1,:),cond_idx_M1(idx_test_M1,:));
        corr_CC(iNet)=mean(r);

        %% prepare for dPCA
        % dPCA Prep
        [~,prepdata,~]=get_prep_exec_after_FR(states,idx_pos,idx_dir,idx_dist,idx_Ncycle);
        %plot_by_dir_and_cycle(endexecdata.scores,endexecdata.ndir,endexecdata.npos,endexecdata.ndist,true)

        figure(figW)
        %% Euclidean distance between trajectories
        [~,~,~,Dist_all_prep(:,:,iNet)]=compare_traj_directions(prepdata.scores(:,1:ndims),prepdata.ndir,prepdata.npos,prepdata.ndist,Npoints);

        if strcmp(ff(iNet).name,'Trained_EMG_Hyp_continuousDrake_14_SLen300v2.mat') ||  strcmp(ff(iNet).name,'Trained_EMG_Hyp_separateDrake_12_Orth_start.mat')  || strcmp(ff(iNet).name,'Trained_M1_Hyp_continuousCousteau_1.mat') || strcmp(ff(iNet).name,'Trained_M1_Hyp_separateCousteau_10_Orth_start.mat')

            states2=states-mean(states);
            states2=states2./repmat(range(states2)+5,size(states,1),1);
            [~,PC_this_example]=pca(states2);
            PC_this_example=PC_this_example-PC_this_example(1,:);
            figure(figW)
            % pos1 = [0.0977 0.5378 0.1212 0.2284];
            % pos2 = [0.2476 0.6666 0.0529 0.0977];
            % pos3 = [0.2448 0.5523 0.0604 0.0956];

            % if ifamily==1
            %     ax=[axes('Position',pos1)...
            %         axes('Position',pos2)...
            %         axes('Position',pos3)];
            % else
            %     pos1(1)=pos1(1)+0.3;
            %     pos2(1)=pos2(1)+0.3;
            %     pos3(1)=pos3(1)+0.3;
            %     ax=[axes('Position',pos1)...
            %         axes('Position',pos2)...
            %         axes('Position',pos3)];
            % end

            subplot(4,4,4+plot_column)
            dt = 10;

            colour_dist=plasma(numel(dists));
            %% plot first and last distances
            for i_dist=1:4:5
                idx=find(idx_dir==2 & idx_pos==1 & idx_dist==dists(i_dist));
                hold on
                plot3(PC_this_example(idx,1),PC_this_example(idx,2),PC_this_example(idx,3),'Color',colour_dist(i_dist,:))
                plot3(PC_this_example(idx(1),1),PC_this_example(idx(1),2),PC_this_example(idx(1),3),'.','MarkerSize',30,'Color',colour_dist(i_dist,:))
                plot3(PC_this_example(idx(100*dt:end-40*dt),1),PC_this_example(idx(100*dt:end-40*dt),2),PC_this_example(idx(100*dt:end-40*dt),3),'Color',colour_dist(i_dist,:),'LineWidth',3)
            end
            view(64.5239,15.9280)

            Inputs_all=[Input,Test_input];
            save(['.\..\..\..\Output_files\Scores_' ff(iNet).name],'PC_this_example','idx_dir','idx_pos','idx_dist','idx_Ncycle','trials_idx','Inputs_all')

            % different hyp
            %     view(-116.1139,40.3981)
            %     ax=gca;
            if plot_supp
                do_plot_pred=1;
            else
                do_plot_pred=0;
            end
        else
            do_plot_pred=0;
        end


        %do_plot_pred=1;

        if  strcmp(ff(iNet).name, 'Trained_M1_Hyp_separateCousteau_10_Orth_start.mat') || strcmp(ff(iNet).name,'Trained_M1_Hyp_continuousCousteau_1.mat')
            figure(figW)
            subplot(4,4,8+plot_column)
            plot_init = 1;
        else
            plot_init = 0;
        end

        [Angle_disc_rhyth(iNet),Init_cond_t(start:start+3,:),Dist2Att]=RNNs_predictions(states,idx_dir,idx_pos,idx_dist,exec2,idx_Ncycle,do_plot_pred,plot_column,plot_init,prep_fig);
        start=start+4;

        if plot_supp && (strcmp(ff(iNet).name, 'Trained_EMG_Hyp_continuousDrake_12_SLen300v2.mat') || strcmp(ff(iNet).name,'Trained_EMG_Hyp_separateDrake_10_Orth_start.mat') || strcmp(ff(iNet).name,'Trained_M1_Hyp_continuousCousteau_1.mat') || strcmp(ff(iNet).name,'Trained_M1_Hyp_separateCousteau_10_Orth_start.mat'))

            figure(dpca_fig)
            do_plot=1;
        else
            do_plot=0;

        end
        [~,Percentage_all_prep(iNet,:)]=dPCA_across_conditions(prepdata.FR-mean(prepdata.FR),prepdata.ndir,prepdata.npos,prepdata.ndist,do_plot,plot_column);

        Dist2Att_all(1:size(Dist2Att,1),iNet,:)=mean(Dist2Att(:,:,:),2); % case N cycle 0.5


    end

    clear Bend
end


SuccesfulNets=sum(R2(:,2)>=0.8);
disp(['Succesfully trained networks = ' num2str(SuccesfulNets) ' of ' num2str(Nnetworks)])


if strcmp(region_name,'M1') && plot_supp && strcmp(animal,'Cousteau')
    %% PLot LDS results
    figure(plot_supp_figs.LDS)
    subplot(2,3,plot_column)
    title(animal)
    epochs2=epochs_all;
    epochs2(epochs_all==0)=nan;
    colour_epoch=[255 255 255;0,90,116;74,148,159;196,196,196; 244,119,127;147,0,58]./255;
    ms=10;
    tstart=-100*ms;
    Ncond=20;
    max_l=round(max(t_end)*ms);
    imagesc([tstart max_l],[1 Ncond],median(epochs2(1:(max_l-tstart)/ms,:,:),3,'omitnan')')
    hold on
    plot([0 0],[0 Ncond],'Color',[0 0 0],'LineWidth',2)
    t_2=reshape([t_end,t_end]',Ncond*2,1)*ms;
    tendexec=reshape([(1:Ncond)-0.5;(1:Ncond)+0.5],Ncond*2,1);
    plot(t_2,tendexec,'Color',[0 0 0],'LineWidth',2)
    colormap(colour_epoch)
    clim([0 5])
    box off
    xlabel('Time to movement onset [ms]')
    ylabel('Condition')
    cbh = colorbar ; %Create Colorbar
    cbh.Ticks = linspace(0.5,4.5,6); %Create 8 ticks from zero to 1
    cbh.TickLabels = {'','only contraction','rotation+ contraction','only rotation','rotation+ expansion','only expansion'} ;
    xlim([tstart max_l+1000])
    figure(figW)
end

corr_CC=[corr_CC,rprep,rexec,rend,corr_CC_control,CC_EMG];


%% define outputs
Results.Angle_rot=Angle_rot;
Results.R2=R2;
Results.Dist_all_prep=Dist_all_prep;
Results.Angle_disc_rhyth=Angle_disc_rhyth;
Results.SuccesfulNets=SuccesfulNets;
Results.corr_CC=corr_CC;
Results.Percentage_all_prep=Percentage_all_prep;
Results.Mean_OutW=Mean_OutW;
Results.Init_cond_t=Init_cond_t;
Results.Dist2Att_all=Dist2Att_all;
end

function [idx_dir,idx_pos,idx_Ncycle,idx_dist,prep,exec,startexec,endexec]=find_idx_conds(trials_idx,idx_conds,idx_current_cycle)

Nconds=max(trials_idx);
ntimepoints=size(trials_idx,1);
idx_dir=zeros(ntimepoints,1);
idx_dist=zeros(ntimepoints,1);
idx_pos=zeros(ntimepoints,1);
idx_Ncycle=nan(ntimepoints,1);
exec=ones(ntimepoints,1);
prep=ones(ntimepoints,1);
endexec=ones(ntimepoints,1);
startexec=zeros(ntimepoints,1);
%NumberCyles=unique(idx_conds(:,1));



for i_cond=1:Nconds
    this_cond=find(trials_idx==i_cond);

    idx_Ncycle(this_cond)=idx_current_cycle(i_cond,1:numel(this_cond));
    idx_pos(this_cond)=idx_conds(i_cond,2);
    idx_dir(this_cond)=idx_conds(i_cond,3);

    exec(this_cond(1:100))=0;
    exec(this_cond(end-40:end))=0;
    % prep idx
    prep(this_cond(101:end))=0;


    idx_mov_end=find(~isnan(idx_current_cycle(i_cond,:)),1,'last');

    if idx_conds(i_cond,1)<1

        end_first_half=idx_mov_end;

    elseif idx_conds(i_cond,1)==1
        firstcycledur=idx_mov_end-100;
        end_first_half=100+round(firstcycledur/2);
    else
        firstcycledur=find(idx_current_cycle(i_cond,100:end)>1,1,"first");
        end_first_half=100+round(firstcycledur/2);

    end
    startexec(this_cond(100:end_first_half))=1;

    % last half cycle
    endexec(this_cond(1:end-100))=0;
    endexec(this_cond(end-40:end))=0;

    %idx_dist(this_cond)=find(idx_conds(i_cond,1)==NumberCyles);
    idx_dist(this_cond)=idx_conds(i_cond,1);
end

end