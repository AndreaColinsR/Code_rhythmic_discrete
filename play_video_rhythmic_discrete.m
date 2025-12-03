function play_video_rhythmic_discrete(region_name,animal,RNN_name_same,RNN_name_diff,i_dir,i_pos)
%% play_video_rhythmic_discrete creates the supplemenary videos 1 and 2
% This function creates a video that shows the neural trajectories of the
% RNNs (RNN_name_same and RNN_name_diff) and the neural recordings from
% the region region_name
%
%% Inputs
%
%  region_name     – Brain region recorded (e.g., 'M1' or 'SMA')
%
%  animal          – Animal identifier (e.g., 'Cousteau' or 'Drake')
%
%  RNN_name_same   – Name of the RNN file modeling same-control movement
%
%  RNN_name_diff   – Name of the RNN file modeling different-control movement
%
%  i_dir           – Movement direction index to visualise (e.g., 1 or 2)
%
%  i_pos           – Movement position index to visualise (e.g., 1 or 2)
% 
%
%% Outputs
% This functions generates frames saved as png files in the folder
% Output_files\Frames, which can be compiled into a video
%
% 03/12/2025
% Andrea Colins Rodriguez


%% Neural recordings
load(['.\Output_files\scores_' animal '_' region_name '.mat'],'scores','idx_dir','idx_pos','idx_dist','exec')
idx_rhythmic=idx_dir==i_dir & idx_pos==i_pos & idx_dist==7;
idx_discrete=idx_dir==i_dir & idx_pos==i_pos & idx_dist==0.5;
scores_rhythmic_M1=scores(idx_rhythmic,:);
scores_discrete_M1=scores(idx_discrete,:);
Nt = size(scores_rhythmic_M1,1);
Nt2 = size(scores_discrete_M1,1);
paramsM1.Nt2 = Nt2;
paramsM1.Nt = Nt;
paramsM1.view=[48,21];

execution=nan(Nt,3);

paramsM1.limx=[min(scores_rhythmic_M1(:,1)) max(scores_rhythmic_M1(:,1))];
paramsM1.limy=[min(scores_rhythmic_M1(:,2)) max(scores_rhythmic_M1(:,2))];
paramsM1.limz=[min(scores_rhythmic_M1(:,3)) max(scores_rhythmic_M1(:,3))];    
paramsM1.title=region_name;

%% RNN Same
[scores_rhythmic_Same, scores_discrete_Same, ~, paramsRNN_same, InputRNN_same] = extract_info_RNN(RNN_name_same,i_dir,i_pos,region_name);
paramsRNN_same.title="Same-control";

%% RNN Different
[scores_rhythmic_diff, scores_discrete_diff, Nt, paramsRNN_diff, InputRNN_diff] = extract_info_RNN(RNN_name_diff,i_dir,i_pos,region_name);
paramsRNN_diff.title="Different-control";
 if strcmp(region_name,'SMA')
     paramsRNN_same.view=[-1,58];
     paramsRNN_diff.view=[136,36];
     paramsM1.view=[-18,4];
     scores_rhythmic_M1(:,2)=-scores_rhythmic_M1(:,2);
     scores_discrete_M1(:,2)=-scores_discrete_M1(:,2);
     paramsM1.limy=-[paramsM1.limy(2) paramsM1.limy(1)];

 end


colour=plasma(5);
params.colour_rhythmic=colour(5,:);
params.colour_discrete=colour(1,:);
fig1=figure;
set(fig1, 'Position',  [0, 100, 1980, 400])


dt=10;
xtime=1:paramsM1.Nt;
xtime_RNN=1:Nt;

counter=1;

%% create frames 
for t=1:Nt
        
    create_frame_case(1,scores_rhythmic_Same,scores_discrete_Same,xtime_RNN,t,InputRNN_same,params,paramsRNN_same)
    create_frame_case(2,scores_rhythmic_diff,scores_discrete_diff,xtime_RNN,t,InputRNN_diff,params,paramsRNN_diff)
    create_frame_case(3,scores_rhythmic_M1,scores_discrete_M1,xtime,t*dt,execution,params,paramsM1)

    %print(['.\Output_files\Frames\' region_name '_' num2str(counter)],'-dpng','-r0')
    counter=counter+1;
    pause(0.01)

    clf
end

close(fig1)
end

function create_frame_case(caseI,scores_rhythmic,scores_discrete,xtime,t,execution_rhythmic,params,paramsM1)
% create_frame_case creates a frame of first 3 PCs of the neural
% trajectories of discrete and rhythmic movement
%% Inputs
%
%  caseI                    – Number of the case that will be plotted (1 = Same-RNN; 2 = Diff-RNN; 3 = Neural recording)
%
%  scores_rhythmic          – Matrix [times, N PC] containing the
%  neural trajectory of the rhythmic case
%
%  scores_discrete          – Matrix [times, N PC] containing the
%  neural trajectory of the discrete case
%
%  xtime                    – Array containing the times of the neural
%  trajectory of the rhythmic case
%
%  t                        – Time of the neural trajectory that will be plotted
%
%  execution_rhythmic       – Inputs of the RNNs that will be plotted
% 
% params                    - Structure containing the parameters of the
% plot such as colours of the trajectories (common for all plots)
%
% paramsM1                  - Structure containing the parameters of the
% plot specific for each case (axes limits, number of times, etc). This
% structure is defined by the function extract_info_RNN
%
% 03/12/2025
% Andrea Colins Rodriguez

%% Define subplots
caseI=caseI*2;
subplots_i =[caseI-[1 0] 6*2+caseI-[1 0] 6*3+caseI-[1 0] 6*3+caseI-[1 0]];
subplot(5,6,subplots_i);

%% Plots
hold on
plot3([paramsM1.limx(1) paramsM1.limx(2)],[0 0],[0 0],'Color',[0.5 0.5 0.5],'Linewidth',2)
plot3([0 0],[paramsM1.limy(1) paramsM1.limy(2)],[0 0],'Color',[0.5 0.5 0.5],'Linewidth',2)
plot3([0 0],[0 0],[paramsM1.limz(1) paramsM1.limz(2)],'Color',[0.5 0.5 0.5],'Linewidth',2)

[~,posx]=max(abs(paramsM1.limx));
text(paramsM1.limx(posx),0,0.1*paramsM1.limx(2),'PC 1','FontSize',14,'Color',[0.5 0.5 0.5]*1)
text(0,paramsM1.limy(2),-0.1*paramsM1.limy(2),'PC 2','FontSize',14,'Color',[0.5 0.5 0.5]*1)
text(-0.1*paramsM1.limz(2),0,paramsM1.limz(2)*1.1,'PC 3','FontSize',14,'Color',[0.5 0.5 0.5]*1)


plot3(scores_rhythmic(1:t,1),scores_rhythmic(1:t,2),scores_rhythmic(1:t,3),'Color',params.colour_rhythmic,'LineWidth',2)

if t<=paramsM1.Nt2

    plot3(scores_discrete(1:t,1),scores_discrete(1:t,2),scores_discrete(1:t,3),'Color',params.colour_discrete,'LineWidth',2)
else
    plot3(scores_discrete(1:paramsM1.Nt2,1),scores_discrete(1:paramsM1.Nt2,2),scores_discrete(1:paramsM1.Nt2,3),'Color',params.colour_discrete,'LineWidth',2)
end

%% Formatting
view(paramsM1.view(1),paramsM1.view(2))
axis off
title(paramsM1.title,'FontSize', 14)

subplot(5,6,6*4+caseI-[1 0]);
hold on
plot(xtime(1:t),execution_rhythmic(1:t,1),'Color','k','LineWidth',2)
plot(xtime(1:t),execution_rhythmic(1:t,2)+2,'Color',params.colour_rhythmic,'LineWidth',2)
plot(xtime(1:t),execution_rhythmic(1:t,3)+2,'Color',params.colour_discrete,'LineWidth',2)
xlim([0 paramsM1.Nt])
ylim([0 3])
axis off

end

function [scores_rhythmic, scores_discrete, Nt, paramsRNN, InputRNN] = extract_info_RNN(RNN_name,i_dir,i_pos,region_name)
% extract_info_RNN extract the informartion of the neural trajectories that
% will be plotted (rhythmic and discrete cases)
%
%% Inputs
%  RNN_name     - Name of the .mat file containing RNN outputs
%
%  i_dir        – Movement direction index to select (e.g., 1 or 2)
%
%  i_pos        – Movement position index to select (e.g., 1 or 2)
%
%  region_name  - Name of brain region ('SMA' or 'M1').
%
%
%% OUTPUTS
%   scores_rhythmic   - Matrix [times, 3 PCs] PC trajectories for rhythmic
%                             trial
%
%   scores_discrete   - Matrix [times, 3 PCs] PC trajectories for discrete
%                             trial
%
%   Nt                - Number of time steps in rhythmic trials.
%
%   paramsRNN         - Structure with plotting and trajectory parameters:
%                               .Nt      — number of rhythmic time steps
%                               .Nt2     — number of discrete time steps
%                               .view    — 3-D view angles for plotting
%                               .limx    — x-axis limits combining both conditions
%                               .limy    — y-axis limits
%                               .limz    — z-axis limits
%
%   InputRNN         - Matrix [times, N inputs] containing extracted input signals for
%                             an example rhythmic and discrete trial.
%
% 
% 03/12/2025
% Andrea Colins Rodriguez


% Load file 
load(['.\Output_files\' RNN_name '.mat'],'PC_this_example','idx_dir','idx_pos','idx_dist','trials_idx','Inputs_all')

% Select trials
idx_rhythmic=find(idx_dir==i_dir & idx_pos==i_pos & idx_dist==7);
idx_discrete=find(idx_dir==i_dir & idx_pos==i_pos & idx_dist==0.5);

% Select trajectories
scores_rhythmic=PC_this_example(idx_rhythmic,:);
scores_discrete=PC_this_example(idx_discrete,:);
Nt = size(scores_rhythmic,1);
Nt2 = size(scores_discrete,1);

% Define parameters for plotting this trajectories of this RNN
paramsRNN.Nt2 =Nt2;
paramsRNN.Nt = Nt;
paramsRNN.view=[-145,38];
paramsRNN.limx=[min([scores_rhythmic(:,1);scores_discrete(:,1)]) max([scores_rhythmic(:,1);scores_discrete(:,1)])];
paramsRNN.limy=[min([scores_rhythmic(:,2);scores_discrete(:,2)]) max([scores_rhythmic(:,2);scores_discrete(:,2)])];
paramsRNN.limz=[min([scores_rhythmic(:,3);scores_discrete(:,3)]) max([scores_rhythmic(:,3);scores_discrete(:,3)])];


%% Define RNN inputs 
itrial=trials_idx(idx_rhythmic(10));
itrialdisc=trials_idx(idx_discrete(10));
if size(Inputs_all,3)>6 && strcmp(region_name,'SMA')

    InputRNN=[squeeze(Inputs_all(:,itrial,[2,6])) Inputs_all(:,itrialdisc,7)];
else
    InputRNN=[squeeze(Inputs_all(:,itrial,[2,3])) Inputs_all(:,itrialdisc,3)];
end

end