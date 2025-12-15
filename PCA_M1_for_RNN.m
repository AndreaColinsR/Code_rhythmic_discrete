function PCA_M1_for_RNN
%% define the inputs and outputs to train an RNN so the output is M1. Hopefully the RNN should look like SMA
animal={'Cousteau','Drake'}; % 'E' or 'F'
timesmov=[-1000 400]; %Ames used 1250:3250 => cycles 2 to 5
region_name='M1';
for i_animal=1:numel(animal)

    %[scores,~,idx_pos,idx_dir,idx_dist,~,~,idx_Ncycle]=extract_trajectories(animal{i_animal},timesmov);
    load(['.\Output_files\scores_' animal{i_animal} '_'  region_name '.mat'],'scores','idx_dir','idx_pos','idx_dist','idx_Ncycle')
    %create_Data_continuous(scores,idx_pos,idx_dir,idx_dist,idx_Ncycle,animal{i_animal},timesmov)
    %create_Data_separate(scores,idx_pos,idx_dir,idx_dist,idx_Ncycle,animal{i_animal},timesmov)

    %% controls 
    %create_Data_separate_control(scores,idx_pos,idx_dir,idx_dist,idx_Ncycle,animal{i_animal},timesmov)
    create_Data_separate2_control(scores,idx_pos,idx_dir,idx_dist,idx_Ncycle,animal{i_animal},timesmov)
    %create_Data_continuous_control(scores,idx_pos,idx_dir,idx_dist,idx_Ncycle,animal{i_animal},timesmov)


end
end



function [scores,explained,idx_pos,idx_dir,idx_dist,baseline,mov_FR,idx_Ncycle]=extract_trajectories(animal,timesmov)
% timesmov = array with timesmov(1) time relative to mov onset timesmov(2)
% time relative to mov end

load(['.\Data_Russo\' animal '_tt.mat'])
if strcmp(animal,'Cousteau')
    P=Pc;
    clear Pc
else
    P=Pd;
    clear Pd
end

Ncycle=P.mask.cycleNum;
Ndistance=P.mask.dist;
condNum=P.mask.condNum;
Ncond=max(condNum);


%%select conditions before doing PCA
%here same direction, same starting positions, different number of cycles
selected_conditions=zeros(Ncond,1);

% find threshold for speed on distance 0.5 condition
speed=zeros(sum(Ndistance==0.5),4);

for i_ex=1:4
    speed(:,i_ex)=sqrt(P.vA(Ndistance==0.5,i_ex*2-1).^2+P.vA(Ndistance==0.5,i_ex*2).^2);
end

Threshold=mean(speed(1500,:));
Exec_half=double(mean(speed,2)>Threshold);
Mov_end=find(diff(Exec_half)==-1);

Exec_half(Exec_half==0)=nan;
Ncycle(Ndistance==0.5)=Exec_half;

for i=1:Ncond
    selected_conditions(i)= any(condNum==i & Ncycle>=1)*i;
end

selected_conditions(selected_conditions<1)=[];
Ncond=numel(selected_conditions);
% plot(Ncycle)
% hold on
% plot(condNum)
% pause


idx_pos=[];
idx_Ncycle=[];
idx_dir=[];
idx_dist=[];
mov_FR=[];
baseline=[];

for i=1:Ncond
    idx=find(condNum==selected_conditions(i));
    baseline=[baseline;P.xA_raw(idx(1:100),:)];

    mov_onset=find(~isnan(Ncycle(idx)),1,'first');


    tstart=mov_onset+timesmov(1);
    tend=find(~isnan(Ncycle(idx)),1,'last')+timesmov(2);

    idx_Ncycle=[idx_Ncycle;Ncycle(idx(tstart:tend))];
    idx_pos=[idx_pos;P.mask.pos(idx(tstart:tend))*2+1];
    idx_dir=[idx_dir;P.mask.dir(idx(tstart:tend))/2+1.5];
    idx_dist=[idx_dist;Ndistance(idx(tstart:tend))];

    mov_FR=[mov_FR;P.xA_raw(idx(tstart:tend),:)];
    
    
%     plot(P.pA(idx(tstart:tend),1))
%     hold on 
%     plot(P.pA(idx(tstart:tend),2))
%     hold off

end


% exclude neurons with zero FR
selected_neurons=sum(mov_FR)>1e-6;
mov_FR=mov_FR(:,selected_neurons);
baseline=mean(baseline(:,selected_neurons),1);

% soft normalise
norm_factor=range(mov_FR)+5;
mov_FR=mov_FR./norm_factor;


%% prep subspace
[coeffs,scores,~,~,explained]=pca(mov_FR);
baseline=(baseline./norm_factor-mean(mov_FR,1))*coeffs;

% recentre activity around region of spontaneous activity
scores=scores-baseline;
end

function create_Data_separate(scores,idx_pos,idx_dir,idx_dist,idx_Ncycle,animal,timesmov)

% this function saves its outputs into a file
scaling=10; %timesteps is 4 ms
Slen=700;
colour_pos=[0 1 1; 1 0 0];
NcondPerCycle=4;


Ncycle=unique(idx_dist);
simLength=round(5000/scaling);% 5000 ms divided by the timesteps
Output=zeros(numel(Ncycle)*2,simLength,4);

%% 7 inputs
% 1 and 2 indicate the start of mov and start pos
% 3 indicates stop signal for rhythmic mov
% 4 and 5 indicate direction
% 6 and 7 indicate if it's rhythmic or discrete

NNcycle=numel(Ncycle);
% for Slen300v2
Input=zeros(NNcycle*NcondPerCycle,simLength,7);
%otherwise
%Input=zeros(NNcycle*NcondPerCycle,simLength,8);
%Input2=zeros(NNcycle*NcondPerCycle,simLength,1); % cotinuous input during execution
idx_dir_trial=zeros(NNcycle*NcondPerCycle,1);
idx_pos_trial=zeros(NNcycle*NcondPerCycle,1);
idx_cycle_trial=zeros(NNcycle*NcondPerCycle,1);


figure
hold on
counter=1;
exec=zeros(NNcycle*NcondPerCycle,simLength);
idx_current_cycle=nan(NNcycle*NcondPerCycle,simLength);

for i_dist=1:numel(Ncycle)
    for direction=1:2
        for i_pos=1:2

            idx=find(idx_dir==direction & idx_pos==i_pos & idx_dist==Ncycle(i_dist));

            subplot(1,2,1)
            hold on
            plot3(scores(idx,1),scores(idx,2),scores(idx,3),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot3(scores(idx(1),1),scores(idx(1),2),scores(idx(1),3),'o','MarkerFaceColor',colour_pos(i_pos,:)./direction,'MarkerEdgeColor',colour_pos(i_pos,:)./direction)

            for i_PC=1:4
                subplot(8,2,i_PC*2)
                hold on
                plot(scores(idx,i_PC),'Color',colour_pos(i_pos,:),'LineWidth',2)
                ylabel(['PC ' num2str(i_PC)])
            end

            % format into tensorflow output [trials,timesteps,Noutputs]
            Ntimes=numel(idx);

            %% Define Outputs as the first 4 component of EMG
            Output(counter,1:round(Ntimes/scaling),:)=interp1(1:Ntimes,scores(idx,1:4),linspace(1,Ntimes,round(Ntimes/scaling)));

            % Set the final steady state
            Output(counter,round(Ntimes/scaling)+1:end,:)=repmat(Output(counter,round(Ntimes/scaling),:),1,simLength-round(Ntimes/scaling),1);
            %% Define inputs

            endt=round((numel(idx)-timesmov(2))/scaling);
            startt=round((-timesmov(1))/scaling);

            %% pos
            Input(counter,startt-Slen/scaling:startt,1)=i_pos-1;
            Input(counter,startt-Slen/scaling:startt,2)=round(-(i_pos-2)/2);
         
            %Input(counter,endt-round(Slen/scaling)-round(100/scaling):endt-round(100/scaling),3)=1;

            %% direction
            Input(counter,startt-Slen/scaling:startt,4)=direction-1;
            Input(counter,startt-Slen/scaling:startt,5)=round(-(direction-2)/2);
            

            %% Defining Decreasing input

            inputL=endt-round(100/scaling)-(startt-round(300/scaling))+1;

            % long input 
            Input(counter,startt-Slen/scaling:startt,6)=1*(Ncycle(i_dist)>=2);
            Input(counter,startt-Slen/scaling:startt,7)=1*(Ncycle(i_dist)<2);

            Input(counter,startt-round(300/scaling):endt-round(100/scaling),6)=linspace(1,0,inputL).*(Ncycle(i_dist)>=2);
            Input(counter,startt-round(300/scaling):endt-round(100/scaling),7)=linspace(1,0,inputL).*(Ncycle(i_dist)<2);

            exec(counter,startt:endt)=1;

            idx_current_cycle(counter,1:round(Ntimes/scaling))=interp1(1:Ntimes,idx_Ncycle(idx),linspace(1,Ntimes,round(Ntimes/scaling)));
            %Input2(counter,round(800/scaling):endt-round(500/scaling),1)=i_pos*2-3;

            % format into tensorflow input [trials,timesteps,Ninputs]
            subplot(8,2,i_PC*2+2)
            hold on
            plot(Output(counter,:,1),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+4)
            hold on
            plot(Input(counter,:,3)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+6)
            hold on
            plot(Input(counter,:,6)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            %plot(Input(counter,:,3)+counter+0.5,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+8)
            hold on
            plot(Input(counter,:,4)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            %plot(Input(counter,:,5)+counter+0.5,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            idx_pos_trial(counter)=i_pos;
            idx_dir_trial(counter)=direction;
            idx_cycle_trial(counter)=Ncycle(i_dist);

            counter=counter+1;
        end
    end
end

xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

%% Initial position input is projected by two random vectors
Bipos=normrnd(0, 0.3,50,1);
Bipos2 = rotate_n_dimensional_vector(Bipos,pi);
Bipos=[Bipos,Bipos2]';

%% direction input is projected by two orthogonal random vectors

Bdir=normrnd(0, 0.3,50,1);
Bdir2 = rotate_n_dimensional_vector(Bdir,pi/2);
Bdir=[Bdir,Bdir2]';

%% task type is projected as two orthogonal random vectors
Btask=normrnd(0, 0.3,50,1);
Btask2 = rotate_n_dimensional_vector(Btask,pi/2);
Btask=[Btask,Btask2]';

Bend=normrnd(0, 0.3,50,1);
Bend2 = rotate_n_dimensional_vector(Bend,pi/2);
Bend=[Bend,Bend2]';

save(['.\Code\Training_RNNs\M1_data\M1_' animal '_separate_long.mat'],'Input','Output','idx_pos_trial','idx_dir_trial','idx_cycle_trial','Bipos','Bdir','Btask','exec','idx_current_cycle')

end

function create_Data_continuous(scores,idx_pos,idx_dir,idx_dist,idx_Ncycle,animal,timesmov)

% this function saves its outputs into a file
scaling=10; %timesteps is 4 ms
Slen=700;
colour_pos=[0 1 1; 1 0 0];
NcondPerCycle=4;

Ncycle=unique(idx_dist);
simLength=round(5000/scaling);% 5000 ms divided by the timesteps
Output=zeros(numel(Ncycle)*2,simLength,4);

%% 5 inputs
% 1 and 2 indicate the start of mov and start pos
% 3 indicates movement end
% 4 and 5 indicate the start of mov (same timing than 1 and 2) and
% direction
NNcycle=numel(Ncycle);
Input=zeros(NNcycle*NcondPerCycle,simLength,5);
Input2=zeros(NNcycle*NcondPerCycle,simLength,1); % cotinuous input during execution
idx_dir_trial=zeros(NNcycle*NcondPerCycle,1);
idx_pos_trial=zeros(NNcycle*NcondPerCycle,1);
idx_cycle_trial=zeros(NNcycle*NcondPerCycle,1);
exec=zeros(NNcycle*NcondPerCycle,simLength);
idx_current_cycle=nan(NNcycle*NcondPerCycle,simLength);


figure
hold on
counter=1;

for i_dist=1:numel(Ncycle)
    for direction=1:2
        for i_pos=1:2

            idx=find(idx_dir==direction & idx_pos==i_pos & idx_dist==Ncycle(i_dist));

            subplot(1,2,1)
            hold on
            plot3(scores(idx,1),scores(idx,2),scores(idx,3),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot3(scores(idx(1),1),scores(idx(1),2),scores(idx(1),3),'o','MarkerFaceColor',colour_pos(i_pos,:)./direction,'MarkerEdgeColor',colour_pos(i_pos,:)./direction)


            for i_PC=1:4
                subplot(8,2,i_PC*2)
                hold on
                plot(scores(idx,i_PC),'Color',colour_pos(i_pos,:),'LineWidth',2)
                ylabel(['PC ' num2str(i_PC)])
            end

            % format into tensorflow output [trials,timesteps,Noutputs]
            Ntimes=numel(idx);

            %% Define Outputs as the first 4 component of EMG
            Output(counter,1:round(Ntimes/scaling),:)=interp1(1:Ntimes,scores(idx,1:4),linspace(1,Ntimes,round(Ntimes/scaling)));

            % Set the final steady state
            Output(counter,round(Ntimes/scaling)+1:end,:)=repmat(Output(counter,round(Ntimes/scaling),:),1,simLength-round(Ntimes/scaling),1);
            %% Define inputs
            endt=round((Ntimes-timesmov(2))/scaling);
            startt=round((-timesmov(1))/scaling);
            % pos
            Input(counter, startt-round(Slen/scaling):startt,1)=i_pos-1;
            Input(counter, startt-round(Slen/scaling):startt,2)=round(-(i_pos-2)/2);


            %% Defining Decreasing input
            inputL=endt-round(100/scaling)-(startt-round(300/scaling))+1;
            
            Input(counter,startt-round(Slen/scaling):startt,3)=1;
            Input(counter,startt-round(300/scaling):endt-round(100/scaling),3)=linspace(1,0,inputL);
            % dir
            Input(counter,startt-round(Slen/scaling):startt,4)=direction-1;
            Input(counter,startt-round(Slen/scaling):startt,5)=round(-(direction-2)/2);


            Input2(counter,round(800/scaling):endt-round(500/scaling),1)=i_pos*2-3;
            exec(counter,startt:endt)=1;
            idx_current_cycle(counter,1:round(Ntimes/scaling))=interp1(1:Ntimes,idx_Ncycle(idx),linspace(1,Ntimes,round(Ntimes/scaling)));
            % format into tensorflow input [trials,timesteps,Ninputs]
            subplot(8,2,i_PC*2+2)
            hold on
            plot(Output(counter,:,1),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot(exec(counter,:))

            subplot(8,2,i_PC*2+4)
            hold on
            plot(Input(counter,:,1)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+6)
            hold on
            %plot(Input(counter,:,2)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot(Input(counter,:,3)+counter+0.5,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+8)
            hold on
            plot(Input(counter,:,4)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot(Input(counter,:,5)+counter+0.5,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            idx_pos_trial(counter)=i_pos;
            idx_dir_trial(counter)=direction;
            idx_cycle_trial(counter)=Ncycle(i_dist);

            counter=counter+1;
        end
    end
end

xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

%% Initial position input is projected by two random vectors
Bipos=rand(50,1)-0.5;
Bipos2 = rotate_n_dimensional_vector(Bipos,pi);
Bipos=[Bipos,Bipos2]';

%% direction input is projected by two orthogonal random vectors

Bdir=rotate_n_dimensional_vector(Bipos2,pi/2);
Bdir2 =rotate_n_dimensional_vector(Bdir,pi);
Bdir=[Bdir,Bdir2]';

save(['.\Code\Training_RNNs\M1_data\M1_' animal '_continuous_long.mat'],'Input','Output','idx_pos_trial','idx_cycle_trial','idx_dir_trial','Bipos','Bdir','exec','idx_current_cycle')
end

function create_Data_continuous_control(scores,idx_pos,idx_dir,idx_dist,idx_Ncycle,animal,timesmov)

% this function saves its outputs into a file
scaling=10; %timesteps is 4 ms
Slen=700;
colour_pos=[0 1 1; 1 0 0];
NcondPerCycle=4;

Ncycle=unique(idx_dist);
simLength=round(5000/scaling);% 5000 ms divided by the timesteps
Output=zeros(numel(Ncycle)*2,simLength,4);

%% 5 inputs
% 1 and 2 indicate the start of mov and start pos
% 3 indicates movement end
% 4 and 5 indicate the start of mov (same timing than 1 and 2) and
% direction
NNcycle=numel(Ncycle);
Input=zeros(NNcycle*NcondPerCycle,simLength,5);
Input2=zeros(NNcycle*NcondPerCycle,simLength,1); % cotinuous input during execution
idx_dir_trial=zeros(NNcycle*NcondPerCycle,1);
idx_pos_trial=zeros(NNcycle*NcondPerCycle,1);
idx_cycle_trial=zeros(NNcycle*NcondPerCycle,1);
exec=zeros(NNcycle*NcondPerCycle,simLength);
idx_current_cycle=nan(NNcycle*NcondPerCycle,simLength);


figure
hold on
counter=1;

for i_dist=1:numel(Ncycle)
    for direction=1:2
        for i_pos=1:2

            idx=find(idx_dir==direction & idx_pos==i_pos & idx_dist==Ncycle(i_dist));

            subplot(1,2,1)
            hold on
            plot3(scores(idx,1),scores(idx,2),scores(idx,3),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot3(scores(idx(1),1),scores(idx(1),2),scores(idx(1),3),'o','MarkerFaceColor',colour_pos(i_pos,:)./direction,'MarkerEdgeColor',colour_pos(i_pos,:)./direction)


            for i_PC=1:4
                subplot(8,2,i_PC*2)
                hold on
                plot(scores(idx,i_PC),'Color',colour_pos(i_pos,:),'LineWidth',2)
                ylabel(['PC ' num2str(i_PC)])
            end

            % format into tensorflow output [trials,timesteps,Noutputs]
            Ntimes=numel(idx);

            %% Define Outputs as the first 4 component of EMG
            Output(counter,1:round(Ntimes/scaling),:)=interp1(1:Ntimes,scores(idx,1:4),linspace(1,Ntimes,round(Ntimes/scaling)));

            % Set the final steady state
            Output(counter,round(Ntimes/scaling)+1:end,:)=repmat(Output(counter,round(Ntimes/scaling),:),1,simLength-round(Ntimes/scaling),1);
            %% Define inputs
            endt=round((Ntimes-timesmov(2))/scaling);
            startt=round((-timesmov(1))/scaling);
            % pos
            Input(counter, startt-round(Slen/scaling):startt,1)=i_pos-1;
            Input(counter, startt-round(Slen/scaling):startt,2)=round(-(i_pos-2)/2);


            %% Defining Decreasing input
            %inputL=endt-round(100/scaling)-(startt-round(300/scaling))+1;
            
            % in control we have same inputs than for M1 network:
            % two inputs
            % 1 input to start the movement (here direction and pos= 4,5)
            % 1 input to stop the movement
            Input(counter,endt-round(300/scaling):endt,3)=1;
            % dir
            Input(counter,startt-round(Slen/scaling):startt,4)=direction-1;
            Input(counter,startt-round(Slen/scaling):startt,5)=round(-(direction-2)/2);


            Input2(counter,round(800/scaling):endt-round(500/scaling),1)=i_pos*2-3;
            exec(counter,startt:endt)=1;
            idx_current_cycle(counter,1:round(Ntimes/scaling))=interp1(1:Ntimes,idx_Ncycle(idx),linspace(1,Ntimes,round(Ntimes/scaling)));
            % format into tensorflow input [trials,timesteps,Ninputs]
            subplot(8,2,i_PC*2+2)
            hold on
            plot(Output(counter,:,1),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot(exec(counter,:))

            subplot(8,2,i_PC*2+4)
            hold on
            plot(Input(counter,:,1)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+6)
            hold on
            %plot(Input(counter,:,2)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot(Input(counter,:,3)+counter+0.5,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+8)
            hold on
            plot(Input(counter,:,4)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot(Input(counter,:,5)+counter+0.5,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            idx_pos_trial(counter)=i_pos;
            idx_dir_trial(counter)=direction;
            idx_cycle_trial(counter)=Ncycle(i_dist);

            counter=counter+1;
        end
    end
end

xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

%% Initial position input is projected by two random vectors
Bipos=rand(50,1)-0.5;
Bipos2 = rotate_n_dimensional_vector(Bipos,pi);
Bipos=[Bipos,Bipos2]';

%% direction input is projected by two orthogonal random vectors

Bdir=rotate_n_dimensional_vector(Bipos2,pi/2);
Bdir2 =rotate_n_dimensional_vector(Bdir,pi);
Bdir=[Bdir,Bdir2]';

save(['.\Code\Training_RNNs\M1_data\M1_' animal '_continuous_long_control.mat'],'Input','Output','idx_pos_trial','idx_cycle_trial','idx_dir_trial','Bipos','Bdir','exec','idx_current_cycle')
end

function create_Data_separate_control(scores,idx_pos,idx_dir,idx_dist,idx_Ncycle,animal,timesmov)

% this function saves its outputs into a file
scaling=10; %timesteps is 4 ms
Slen=700;
colour_pos=[0 1 1; 1 0 0];
NcondPerCycle=4;


Ncycle=unique(idx_dist);
simLength=round(5000/scaling);% 5000 ms divided by the timesteps
Output=zeros(numel(Ncycle)*2,simLength,4);

%% 7 inputs
% 1 and 2 indicate the start of mov and start pos
% 3 indicates movement end for rhythmic mov
% 4 and 5 indicate the start of mov (same timing than 1 and 2) and
% direction
% 6 and 7 indicate if it's rhythmic or discrete

NNcycle=numel(Ncycle);
% for Slen300v2
Input=zeros(NNcycle*NcondPerCycle,simLength,7);
%otherwise
%Input=zeros(NNcycle*NcondPerCycle,simLength,8);
%Input2=zeros(NNcycle*NcondPerCycle,simLength,1); % cotinuous input during execution
idx_dir_trial=zeros(NNcycle*NcondPerCycle,1);
idx_pos_trial=zeros(NNcycle*NcondPerCycle,1);
idx_cycle_trial=zeros(NNcycle*NcondPerCycle,1);


figure
hold on
counter=1;
exec=zeros(NNcycle*NcondPerCycle,simLength);
idx_current_cycle=nan(NNcycle*NcondPerCycle,simLength);

for i_dist=1:numel(Ncycle)
    for direction=1:2
        for i_pos=1:2

            idx=find(idx_dir==direction & idx_pos==i_pos & idx_dist==Ncycle(i_dist));

            subplot(1,2,1)
            hold on
            plot3(scores(idx,1),scores(idx,2),scores(idx,3),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot3(scores(idx(1),1),scores(idx(1),2),scores(idx(1),3),'o','MarkerFaceColor',colour_pos(i_pos,:)./direction,'MarkerEdgeColor',colour_pos(i_pos,:)./direction)

            for i_PC=1:4
                subplot(8,2,i_PC*2)
                hold on
                plot(scores(idx,i_PC),'Color',colour_pos(i_pos,:),'LineWidth',2)
                ylabel(['PC ' num2str(i_PC)])
            end

            % format into tensorflow output [trials,timesteps,Noutputs]
            Ntimes=numel(idx);

            %% Define Outputs as the first 4 component of EMG
            Output(counter,1:round(Ntimes/scaling),:)=interp1(1:Ntimes,scores(idx,1:4),linspace(1,Ntimes,round(Ntimes/scaling)));

            % Set the final steady state
            Output(counter,round(Ntimes/scaling)+1:end,:)=repmat(Output(counter,round(Ntimes/scaling),:),1,simLength-round(Ntimes/scaling),1);
            %% Define inputs

            endt=round((numel(idx)-timesmov(2))/scaling);
            startt=round((-timesmov(1))/scaling);

            %% pos
            Input(counter,startt-Slen/scaling:startt,1)=i_pos-1;
            Input(counter,startt-Slen/scaling:startt,2)=round(-(i_pos-2)/2);
         
            %Input(counter,endt-round(Slen/scaling)-round(100/scaling):endt-round(100/scaling),3)=1;

            %% direction
            Input(counter,startt-Slen/scaling:startt,4)=direction-1;
            Input(counter,startt-Slen/scaling:startt,5)=round(-(direction-2)/2);
            

            %% Defining Decreasing input

            %inputL=endt-round(100/scaling)-(startt-round(300/scaling))+1;

            % long input 
            % prep
            Input(counter,startt-Slen/scaling:startt,6)=1*(Ncycle(i_dist)>=2);
            Input(counter,startt-Slen/scaling:startt,7)=1*(Ncycle(i_dist)<2);
            
            %Input(counter,endt-round(300/scaling):endt,6)=(Ncycle(i_dist)>=2);
            %Input(counter,endt-round(300/scaling):endt,7)=(Ncycle(i_dist)<2);

            exec(counter,startt:endt)=1;

            idx_current_cycle(counter,1:round(Ntimes/scaling))=interp1(1:Ntimes,idx_Ncycle(idx),linspace(1,Ntimes,round(Ntimes/scaling)));
            %Input2(counter,round(800/scaling):endt-round(500/scaling),1)=i_pos*2-3;

            % format into tensorflow input [trials,timesteps,Ninputs]
            subplot(8,2,i_PC*2+2)
            hold on
            plot(Output(counter,:,1),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+4)
            hold on
            plot(Input(counter,:,3)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+6)
            hold on
            plot(Input(counter,:,6)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            %plot(Input(counter,:,3)+counter+0.5,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+8)
            hold on
            plot(Input(counter,:,4)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            %plot(Input(counter,:,5)+counter+0.5,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            idx_pos_trial(counter)=i_pos;
            idx_dir_trial(counter)=direction;
            idx_cycle_trial(counter)=Ncycle(i_dist);

            counter=counter+1;
        end
    end
end

xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

%% Initial position input is projected by two random vectors
Bipos=normrnd(0, 0.3,50,1);
Bipos2 = rotate_n_dimensional_vector(Bipos,pi);
Bipos=[Bipos,Bipos2]';

%% direction input is projected by two orthogonal random vectors

Bdir=normrnd(0, 0.3,50,1);
Bdir2 = rotate_n_dimensional_vector(Bdir,pi/2);
Bdir=[Bdir,Bdir2]';

%% task type is projected as two orthogonal random vectors
Btask=normrnd(0, 0.3,50,1);
Btask2 = rotate_n_dimensional_vector(Btask,pi/2);
Btask=[Btask,Btask2]';

Bend=normrnd(0, 0.3,50,1);
Bend2 = rotate_n_dimensional_vector(Bend,pi/2);
Bend=[Bend,Bend2]';

save(['.\Code\Training_RNNs\M1_data\M1_' animal '_separate_long_control.mat'],'Input','Output','idx_pos_trial','idx_dir_trial','idx_cycle_trial','Bipos','Bdir','Btask','exec','idx_current_cycle')

end

function create_Data_separate2_control(scores,idx_pos,idx_dir,idx_dist,idx_Ncycle,animal,timesmov)

% this function saves its outputs into a file
scaling=10; %timesteps is 4 ms
Slen=700;
colour_pos=[0 1 1; 1 0 0];
NcondPerCycle=4;


Ncycle=unique(idx_dist);
simLength=round(5000/scaling);% 5000 ms divided by the timesteps
Output=zeros(numel(Ncycle)*2,simLength,4);

%% 7 inputs
% 1 and 2 indicate the start of mov and start pos
% 3 indicates movement end for rhythmic mov
% 4 and 5 indicate the start of mov (same timing than 1 and 2) and
% direction
% 6 and 7 indicate if it's rhythmic or discrete

NNcycle=numel(Ncycle);
% for Slen300v2
Input=zeros(NNcycle*NcondPerCycle,simLength,7);
%otherwise
%Input=zeros(NNcycle*NcondPerCycle,simLength,8);
%Input2=zeros(NNcycle*NcondPerCycle,simLength,1); % cotinuous input during execution
idx_dir_trial=zeros(NNcycle*NcondPerCycle,1);
idx_pos_trial=zeros(NNcycle*NcondPerCycle,1);
idx_cycle_trial=zeros(NNcycle*NcondPerCycle,1);


figure
hold on
counter=1;
exec=zeros(NNcycle*NcondPerCycle,simLength);
idx_current_cycle=nan(NNcycle*NcondPerCycle,simLength);

for i_dist=1:numel(Ncycle)
    for direction=1:2
        for i_pos=1:2

            idx=find(idx_dir==direction & idx_pos==i_pos & idx_dist==Ncycle(i_dist));

            subplot(1,2,1)
            hold on
            plot3(scores(idx,1),scores(idx,2),scores(idx,3),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            plot3(scores(idx(1),1),scores(idx(1),2),scores(idx(1),3),'o','MarkerFaceColor',colour_pos(i_pos,:)./direction,'MarkerEdgeColor',colour_pos(i_pos,:)./direction)

            for i_PC=1:4
                subplot(8,2,i_PC*2)
                hold on
                plot(scores(idx,i_PC),'Color',colour_pos(i_pos,:),'LineWidth',2)
                ylabel(['PC ' num2str(i_PC)])
            end

            % format into tensorflow output [trials,timesteps,Noutputs]
            Ntimes=numel(idx);

            %% Define Outputs as the first 4 component of EMG
            Output(counter,1:round(Ntimes/scaling),:)=interp1(1:Ntimes,scores(idx,1:4),linspace(1,Ntimes,round(Ntimes/scaling)));

            % Set the final steady state
            Output(counter,round(Ntimes/scaling)+1:end,:)=repmat(Output(counter,round(Ntimes/scaling),:),1,simLength-round(Ntimes/scaling),1);
            %% Define inputs

            endt=round((numel(idx)-timesmov(2))/scaling);
            startt=round((-timesmov(1))/scaling);

            %% pos
            Input(counter,startt-Slen/scaling:startt,1)=i_pos-1;
            Input(counter,startt-Slen/scaling:startt,2)=round(-(i_pos-2)/2);

            Input(counter,endt-round(400/scaling):endt-round(100/scaling),3)=1;
            %Input(counter,endt-round(Slen/scaling)-round(100/scaling):endt-round(100/scaling),3)=1;

            %% direction
            Input(counter,startt-Slen/scaling:startt,4)=direction-1;
            Input(counter,startt-Slen/scaling:startt,5)=round(-(direction-2)/2);
            

            %% Defining Decreasing input
            %inputL=endt-round(100/scaling)-(startt-round(300/scaling))+1;

            % long input 
            % prep
            Input(counter,startt-Slen/scaling:startt,6)=1*(Ncycle(i_dist)>=2);
            Input(counter,startt-Slen/scaling:startt,7)=1*(Ncycle(i_dist)<2);
            
            %Input(counter,endt-round(300/scaling):endt,6)=(Ncycle(i_dist)>=2);
            %Input(counter,endt-round(300/scaling):endt,7)=(Ncycle(i_dist)<2);

            exec(counter,startt:endt)=1;

            idx_current_cycle(counter,1:round(Ntimes/scaling))=interp1(1:Ntimes,idx_Ncycle(idx),linspace(1,Ntimes,round(Ntimes/scaling)));
            %Input2(counter,round(800/scaling):endt-round(500/scaling),1)=i_pos*2-3;

            % format into tensorflow input [trials,timesteps,Ninputs]
            subplot(8,2,i_PC*2+2)
            hold on
            plot(Output(counter,:,1),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+4)
            hold on
            plot(Input(counter,:,3)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+6)
            hold on
            plot(Input(counter,:,6)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            %plot(Input(counter,:,3)+counter+0.5,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+8)
            hold on
            plot(Input(counter,:,4)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
            %plot(Input(counter,:,5)+counter+0.5,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            idx_pos_trial(counter)=i_pos;
            idx_dir_trial(counter)=direction;
            idx_cycle_trial(counter)=Ncycle(i_dist);

            counter=counter+1;
        end
    end
end

xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

%% Initial position input is projected by two random vectors
Bipos=normrnd(0, 0.3,50,1);
Bipos2 = rotate_n_dimensional_vector(Bipos,pi);
Bipos=[Bipos,Bipos2]';

%% direction input is projected by two orthogonal random vectors

Bdir=normrnd(0, 0.3,50,1);
Bdir2 = rotate_n_dimensional_vector(Bdir,pi/2);
Bdir=[Bdir,Bdir2]';

%% task type is projected as two orthogonal random vectors
Btask=normrnd(0, 0.3,50,1);
Btask2 = rotate_n_dimensional_vector(Btask,pi/2);
Btask=[Btask,Btask2]';

Bend=normrnd(0, 0.3,50,1);
Bend2 = rotate_n_dimensional_vector(Bend,pi/2);
Bend=[Bend,Bend2]';

%save(['.\Code\Training_RNNs\M1_data\M1_' animal '_separate2_long_control.mat'],'Input','Output','idx_pos_trial','idx_dir_trial','idx_cycle_trial','Bipos','Bdir','Btask','exec','idx_current_cycle')

end