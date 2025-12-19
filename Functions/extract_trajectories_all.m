function [scores,explained,idx_pos,idx_dir,idx_dist,baseline,mov_FR,idx_Ncycle,exec]=extract_trajectories_all(animal,region,timesmov,dataset_path)
%% extract_trajectories_all performs PCA on the neural activity of a recording and define neural trajectories for all behavioural conditions
%
%% Inputs
%
% animal: Name of the animal to be analysed. e.g 'Drake' or 'Cousteau'
%
% region: Name of the region to be analysed. e.g 'M1', 'SMA' or 'EMG'
%
% timesmov: array indicating the times to select behavioural and neural
% data. Times selected will start at timesmov(1) before movement onset up
% to timesmov(2) after movement offset
%
% dataset_path: folder path where recordings are stored
%
%% Outputs
%
% scores: Matrix [N_times x 10] of neural activity projected onto the first 10 Principal Components (PCs)
%
% explained: Array [N_neurons] indicating the variance explained by each
% PC
%
% idx_pos: array of N_Times elementes indicating the starting position of the trial: bottom start (1), top start (2)
%
% idx_dir: array of N_Times elementes denoting if in a given timebin corresponded to a trial when the animal cycled forward (2) or backward (1)
%
% idx_dist: array of N_Times elementes indicating the pedalling distance of the trial (0.5, 1, 2, 4, 7)
%
% baseline: Projection of the average neural activity between 1000-900 ms before movement onset
% on the neural subspace
%
% mov_FR: Matrix [N_times x N_neurons] Average firing rate of each neuron for
% each trial type. FR published was already filtered with a gaussian filter
% (sigma = 25).
%
% idx_Ncycle: array of N_Times elementes indicating the number of the current cycle performed (1 to 7, Nan if not executing)
%
% exec: array of N_Times elementes denoting if in a given timebin the animal was executing a movement (1) or not (0)
%
%
% 04/09/2025
% Andrea Colins Rodriguez

load([dataset_path '\' animal '_tt.mat'])

if strcmp(animal,'Cousteau')
    P=Pc;
    clear Pc
else
    P=Pd;
    clear Pd
end


%% M1 (‘.xA_raw’)
%% SMA (‘.uA_raw’)
%% EMG (‘.zA_raw’)

if strcmp(region,'M1')
    P.xA_raw=P.xA_raw;
elseif strcmp(region,'SMA')
    P.xA_raw=P.uA_raw;
elseif strcmp(region,'EMG')
    P.xA_raw=P.zA_raw;
else
    disp('Region not recorgnised (SMA, M1 or EMG)')
    keyboard
end


Ncycle=P.mask.cycleNum;
Ndistance=P.mask.dist;
condNum=P.mask.condNum;
Ncond=max(condNum);

% Execution times and current number of cycles were not defined in the
% original dataset for Ndist = 0.5. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find a threshold for speed in distance = 0.5 condition
speed=zeros(sum(Ndistance==0.5),4);

% compute speed profile for the 4 examples of speed in the data
for i_ex=1:4
    speed(:,i_ex)=sqrt(P.vA(Ndistance==0.5,i_ex*2-1).^2+P.vA(Ndistance==0.5,i_ex*2).^2);
end

% define mean (across trials) of speed at mov onset for one direction and
% pos. Movement offset is defined as the moment when the speed is lower than the threshold. 
Threshold=mean(speed(1500,:));
Exec_half=double(mean(speed,2)>Threshold);

Exec_half(Exec_half==0)=nan;
Ncycle(Ndistance==0.5)=Exec_half;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_pos=[];
idx_Ncycle=[];
idx_dir=[];
idx_dist=[];
mov_FR=[];
exec=[];

baseline=nan(100*Ncond,size(P.xA_raw,2));

for i=1:Ncond
    idx=find(condNum==i);
    baseline(100*(i-1)+1:100*i,:)=P.xA_raw(idx(1:100),:);

    mov_onset=find(~isnan(Ncycle(idx)),1,'first');

    tstart=mov_onset+timesmov(1);
    tend=find(~isnan(Ncycle(idx)),1,'last')+timesmov(2);

    idx_Ncycle=[idx_Ncycle;Ncycle(idx(tstart:tend))];
    idx_pos=[idx_pos;P.mask.pos(idx(tstart:tend))*2+1];
    idx_dir=[idx_dir;P.mask.dir(idx(tstart:tend))/2+1.5];
    idx_dist=[idx_dist;Ndistance(idx(tstart:tend))];
    tmp_exec=zeros(tend-tstart+1,1);
    tmp_exec(-timesmov(1):end-timesmov(2))=1;
    exec=[exec;tmp_exec];

    mov_FR=[mov_FR;P.xA_raw(idx(tstart:tend),:)];

end


% exclude neurons with zero FR
selected_neurons=sum(mov_FR)>1e-6;
mov_FR=mov_FR(:,selected_neurons);
baseline=mean(baseline(:,selected_neurons),1);

% soft normalise
norm_factor=range(mov_FR)+5;
mov_FR=mov_FR./norm_factor;


%% pca subspace
[coeffs,scores,~,~,explained]=pca(mov_FR);

% project baseline onto subspace
baseline=(baseline./norm_factor-mean(mov_FR,1))*coeffs;

end