function PCA_EMG_Hyp_separate
animal={'Cousteau','Drake'}; % 'E' or 'F'
timesmov=[-1000 400]; %Ames used 1250:3250 => cycles 2 to 5
%% Ignoring the last cycle improves the performance for a lot!!
%% Be careful in this part before applying
ms=1000;

colour_dir=[[0 0 1];[0 1 0]];
do_extra_plot=0;
all_eigs=[];
mov_dur2=[];

nrec=0;

for i_animal=1:numel(animal)
    Max_eigen_all=[];
    idx_dir_new_all=[];
    idx_pos_new_all=[];
    idx_cycles_all=[];
    mov_exec_all=[];

    nrec=nrec+1;
    [scores,explained,idx_pos,idx_dir,idx_dist,baseline,FR,idx_Ncycle]=extract_trajectories(animal{i_animal},timesmov);



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
%color_dist=plasma(max(Ncycle)+1);
%% load neurons from M1
%total_matrix=Pd.zA;
%px=P.pA(:,1);
%py=P.pA(:,2);
%Cycle_per_trial=[4 7];
%color_dir=[0 0 0; 1 0 0]; %black=forward and red=backward
%color_pos=[0 0 1; 0 1 0]; %blue=0 and gren=0.5

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
    baseline=[baseline;P.zA_raw(idx(1:100),:)];

    mov_onset=find(~isnan(Ncycle(idx)),1,'first');


    tstart=mov_onset+timesmov(1);
    tend=find(~isnan(Ncycle(idx)),1,'last')+timesmov(2);

    idx_Ncycle=[idx_Ncycle;Ncycle(idx(tstart:tend))];
    idx_pos=[idx_pos;P.mask.pos(idx(tstart:tend))*2+1];
    idx_dir=[idx_dir;P.mask.dir(idx(tstart:tend))/2+1.5];
    idx_dist=[idx_dist;Ndistance(idx(tstart:tend))];

    mov_FR=[mov_FR;P.zA_raw(idx(tstart:tend),:)];

    

    subplot(3,1,1)
    plot(mean(P.xA_raw(idx(tstart:tend),:),2),'k')
    hold on
    title('Cortical activity')
    box off
   
    subplot(3,1,2)
    plot(mean(P.zA_raw(idx(tstart:tend),:),2),'-k')
    hold on 
    title('EMG')
    box off

    subplot(3,1,3)
    plot(sqrt(sum(P.vA(idx(tstart:tend),1:2).^2,2)),'k')
    hold on
    title('Hand Speed')
    box off
    
    exec2=zeros(size(idx(tstart:tend)));
    exec2(-timesmov(1):end-timesmov(2))=1;
    plot(exec2/40)
    box off
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
Ncycle=unique(Ndistance);
scaling=10; %timesteps is 4 ms
simLength=round(5000/scaling);% 5000 ms divided by the timesteps
Output=zeros(numel(Ncycle)*2,simLength,4);

%% 7 inputs
% 1 and 2 indicate the start of mov and start pos
% 3 indicates movement end for rhythmic mov
% 4 and 5 indicate the start of mov (same timing than 1 and 2) and
% direction
% 6 and 7 indicate if it's rhythmic or discrete

NcondPerCycle=4;
NNcycle=numel(Ncycle);
% for Slen300v2
Input=zeros(NNcycle*NcondPerCycle,simLength,7);
%otherwise
%Input=zeros(NNcycle*NcondPerCycle,simLength,8);
Input2=zeros(NNcycle*NcondPerCycle,simLength,1); % cotinuous input during execution
idx_dir_trial=zeros(NNcycle*NcondPerCycle,1);
idx_pos_trial=zeros(NNcycle*NcondPerCycle,1);
idx_cycle_trial=zeros(NNcycle*NcondPerCycle,1);
Slen=300;
colour_pos=[0 1 1; 1 0 0];

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
            plot3(scores(idx,1),scores(idx,2),scores(idx,3),'Color',colour_pos(i_pos,:),'LineWidth',2)
            plot3(scores(idx(1),1),scores(idx(1),2),scores(idx(1),3),'o','MarkerFaceColor',colour_pos(i_pos,:)./direction,'MarkerEdgeColor',colour_pos(i_pos,:)./direction)
            hold on
            plot3(baseline(:,1),baseline(:,2),baseline(:,3),'ok')

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

            % pos
            Input(counter,startt-Slen/scaling:startt,1)=i_pos-1;
            Input(counter,startt-Slen/scaling:startt,2)=round(-(i_pos-2)/2);
            

            % for Slen300v2
            Input(counter,endt-round(Slen/scaling)-round(100/scaling):endt-round(100/scaling),3)=1;
            Input(counter,startt-Slen/scaling:startt,4)=direction-1;
            Input(counter,startt-Slen/scaling:startt,5)=round(-(direction-2)/2);
            Input(counter,startt-Slen/scaling:startt,6)=Ncycle(i_dist)>=2;
            Input(counter,startt-Slen/scaling:startt,7)=Ncycle(i_dist)<2;


            % dir
            %Input(counter,startt-Slen/scaling:startt,5)=direction-1;
            %Input(counter,startt-Slen/scaling:startt,6)=round(-(direction-2)/2);


            % task type, this doesn't need to be in the same timing, maybe
            % something brief before all other pulses

            % same start
            %Input(counter,startt-Slen/scaling:startt,7)=1;
            %otherwise
            %Input(counter,startt-Slen/scaling:startt,7)=Ncycle(i_dist)>=2;
            %Input(counter,startt-Slen/scaling:startt,8)=Ncycle(i_dist)<2;

            % end of mov
            %Input(counter,endt-round(Slen/scaling)-round(100/scaling):endt-round(100/scaling),3)=Ncycle(i_dist)>=2;
            %Input(counter,endt-round(Slen/scaling)-round(100/scaling):endt-round(100/scaling),4)=Ncycle(i_dist)<2;
            


            exec(counter,startt:endt)=1;
            idx_current_cycle(counter,1:round(Ntimes/scaling))=interp1(1:Ntimes,idx_Ncycle(idx),linspace(1,Ntimes,round(Ntimes/scaling)));
            Input2(counter,round(800/scaling):endt-round(500/scaling),1)=i_pos*2-3;

            % format into tensorflow input [trials,timesteps,Ninputs]
            subplot(8,2,i_PC*2+2)
            hold on
            plot(Output(counter,:,1),'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+4)
            hold on
            plot(Input(counter,:,7)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)

            subplot(8,2,i_PC*2+6)
            hold on
            plot(Input(counter,:,7)+counter,'Color',colour_pos(i_pos,:)./direction,'LineWidth',2)
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

%subspace(Btask(:,1),Btask(:,2))*180/pi


Bend=normrnd(0, 0.3,50,1);
Bend2 = rotate_n_dimensional_vector(Bend,pi/2);
Bend=[Bend,Bend2]';

subspace(Bend(:,1),Bend(:,2))*180/pi


%save(['.\Code\Training_RNNs\EMG_data\EMG_' animal '_dir_pos_separate_orth.mat'],'Input','Output','idx_pos_trial','idx_dir_trial','idx_cycle_trial','Bipos','Bdir','Btask','exec','Bend','idx_current_cycle')
%save(['.\EMG_' animal '_dir_pos_separate.mat'],'Input','Output','idx_pos_trial','idx_dir_trial','idx_cycle_trial','Bipos','Bdir','Btask','exec','Bend','idx_current_cycle')


save(['.\Code\Training_RNNs\EMG_data\EMG_' animal '_dir_pos_separate_SLen300v2.mat'],'Input','Output','idx_pos_trial','idx_dir_trial','idx_cycle_trial','Bipos','Bdir','Btask','exec','idx_current_cycle')
%Input=Input2;
%save(['.\Code\Training_RNNs\EMG_data\EMG_' animal '_dir_pos_continuous.mat'],'Input','Output','idx_pos_trial','idx_cycle_trial','Bipos')
end