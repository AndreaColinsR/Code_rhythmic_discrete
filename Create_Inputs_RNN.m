function Create_Inputs_RNN

animal={'Cousteau','Drake'}; % 'E' or 'F'
timesmov=[-1000 400];

%% For M1 models

region_name = 'M1';
output_region = 'EMG';

for i_animal=1:numel(animal)
    load(['.\Output_files\scores_' animal{i_animal} '_' output_region '.mat'],'scores', 'idx_dir','idx_Ncycle','idx_pos','idx_dist')

    Create_Inputs_same_M1(animal{i_animal},timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist);
    Create_Inputs_different_M1(animal{i_animal},timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist);

end

%% For SMA models
region_name = 'SMA';
output_region = 'M1';

for i_animal=1:numel(animal)
    load(['.\Output_files\scores_' animal{i_animal} '_' output_region '.mat'],'scores', 'idx_dir','idx_Ncycle','idx_pos','idx_dist')

    Create_Inputs_same_SMA(animal{i_animal},timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist);
    Create_Inputs_different_SMA(animal{i_animal},timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist);

end

end

function Create_Inputs_same_M1(animal,timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist)
%% Create_Inputs_same creates and saves the inputs to the RNNs that model M1 activity same control hypothesis
%% This model has 3 inputs: Start position, movement end and movement direction. Start position and movement direction are one-hot coded so that there are 5 inputs for each trial

%% 5 inputs
% 1 and 2 indicate the start of mov and start pos
% 3 indicates movement end
% 4 and 5 indicate the start of mov (same timing than 1 and 2) and
% direction

%% Parameters
Ncycle=unique(idx_dist); %[0.5, 1,2, 4,7]
NNcycle=numel(Ncycle); % 5 number of cycles
Ninputs = 5;
scaling=10; %timesteps is 10 ms
timesteps=round(5000/scaling);% 5000 ms divided by the timesteps
Ndim_output = 4; % number of dimensions of EMG data
NcondPerCycle = 4; % number of condition for each number of cycles (2 dir * 2 pos)
Ntrials = NNcycle*NcondPerCycle;
Length_stimulus = 300;% ms


%% Initialising structures
% format as in tensorflow output [trials,timesteps,Noutputs]
Output=zeros(Ntrials,timesteps,Ndim_output);

Input=zeros(Ntrials,timesteps,Ninputs);
idx_dir_trial=zeros(Ntrials,1);
idx_pos_trial=zeros(Ntrials,1);
idx_cycle_trial=zeros(Ntrials,1);
exec=zeros(Ntrials,timesteps);
idx_current_cycle=nan(Ntrials,timesteps);

counter = 1;

for i_dist = 1:NNcycle
    for direction=1:2
        for i_pos=1:2

            idx_this_condition=find(idx_dir==direction & idx_pos==i_pos & idx_dist==Ncycle(i_dist));

            % format into tensorflow output [trials,timesteps,Noutputs]
            Ntimes=numel(idx_this_condition);

            %% Define Outputs as the first 4 component of EMG
            Output(counter,1:round(Ntimes/scaling),:)=interp1(1:Ntimes,scores(idx_this_condition,1:Ndim_output),linspace(1,Ntimes,round(Ntimes/scaling)));

            % Set the final steady state
            Output(counter,round(Ntimes/scaling)+1:end,:)=repmat(Output(counter,round(Ntimes/scaling),:),1,timesteps-round(Ntimes/scaling),1);
            
            %% Define inputs
            endt = round((Ntimes-timesmov(2))/scaling); % time of movement offset
            startt = round((-timesmov(1))/scaling); % time of movement onset 

            % Inputs 1 and 2: starting position
            Input(counter, startt-round(Length_stimulus/scaling):startt,1)=i_pos-1;
            Input(counter, startt-round(Length_stimulus/scaling):startt,2)=round(-(i_pos-2)/2);

            % Input 3: end of mov (same for all conditions)
            Input(counter,endt-round(Length_stimulus/scaling)-round(100/scaling):endt-round(100/scaling),3)=1;


            % Inputs 4 and 5: movement direction
            Input(counter,startt-round(Length_stimulus/scaling):startt,4)=direction-1;
            Input(counter,startt-round(Length_stimulus/scaling):startt,5)=round(-(direction-2)/2);
            
            % Execution times (Necesary for post processing of RNNs)
            exec(counter,startt:endt) = 1;
            % current cycle being performed (Necesary for post processing of RNNs)
            idx_current_cycle(counter,1:round(Ntimes/scaling)) = interp1(1:Ntimes,idx_Ncycle(idx_this_condition),linspace(1,Ntimes,round(Ntimes/scaling)));

            % identification of the conditions 
            idx_pos_trial(counter) = i_pos;
            idx_dir_trial(counter) = direction;
            idx_cycle_trial(counter) = Ncycle(i_dist);

            counter=counter+1;
        end
    end
end

save(['.\Output_files\RNNs_Inputs\M1_' animal '_same.mat'],'Input','Output','idx_pos_trial','idx_cycle_trial','idx_dir_trial','exec','idx_current_cycle')
end

function Create_Inputs_different_M1(animal,timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist)
%% Create_Inputs_same creates and saves the inputs to the RNNs that model M1 activity same control hypothesis
%% This model has 3 inputs: Start position, movement end and movement direction. Start position and movement direction are one-hot coded so that there are 5 inputs for each trial

%% 5 inputs
% 1 and 2 indicate the start of mov and start pos
% 3 indicates movement end
% 4 and 5 indicate the start of mov (same timing than 1 and 2) and
% direction

%% Parameters
Ncycle=unique(idx_dist); %[0.5, 1,2, 4,7]
NNcycle=numel(Ncycle); % 5 number of cycles
Ninputs = 7;
scaling=10; %timesteps is 10 ms
timesteps=round(5000/scaling);% 5000 ms divided by the timesteps
Ndim_output = 4; % number of dimensions of EMG data
NcondPerCycle = 4; % number of condition for each number of cycles (2 dir * 2 pos)
Ntrials = NNcycle*NcondPerCycle;
Length_stimulus = 300;% ms


%% Initialising structures
% format as in tensorflow output [trials,timesteps,Noutputs]
Output=zeros(Ntrials,timesteps,Ndim_output);

Input=zeros(Ntrials,timesteps,Ninputs);
idx_dir_trial=zeros(Ntrials,1);
idx_pos_trial=zeros(Ntrials,1);
idx_cycle_trial=zeros(Ntrials,1);
exec=zeros(Ntrials,timesteps);
idx_current_cycle=nan(Ntrials,timesteps);

counter = 1;

for i_dist = 1:NNcycle
    for direction=1:2
        for i_pos=1:2

            idx_this_condition=find(idx_dir==direction & idx_pos==i_pos & idx_dist==Ncycle(i_dist));

            % format into tensorflow output [trials,timesteps,Noutputs]
            Ntimes=numel(idx_this_condition);

            %% Define Outputs as the first 4 component of EMG
            Output(counter,1:round(Ntimes/scaling),:)=interp1(1:Ntimes,scores(idx_this_condition,1:Ndim_output),linspace(1,Ntimes,round(Ntimes/scaling)));

            % Set the final steady state
            Output(counter,round(Ntimes/scaling)+1:end,:)=repmat(Output(counter,round(Ntimes/scaling),:),1,timesteps-round(Ntimes/scaling),1);
            
            %% Define inputs
            endt = round((Ntimes-timesmov(2))/scaling); % time of movement offset
            startt = round((-timesmov(1))/scaling); % time of movement onset 

            % Inputs 1 and 2: starting position
            Input(counter, startt-round(Length_stimulus/scaling):startt,1)=i_pos-1;
            Input(counter, startt-round(Length_stimulus/scaling):startt,2)=round(-(i_pos-2)/2);

            % Input 3: end of mov (same for all conditions)
            Input(counter,endt-round(Length_stimulus/scaling):endt,3)=1;
            

            % Inputs 4 and 5: movement direction
            Input(counter,startt-round(Length_stimulus/scaling):startt,4)=direction-1;
            Input(counter,startt-round(Length_stimulus/scaling):startt,5)=round(-(direction-2)/2);

            % Inputs 6 and 7: indicating the type of movement (key difference with same control model)
            Input(counter,startt-Length_stimulus/scaling:startt,6)=Ncycle(i_dist)>=2;
            Input(counter,startt-Length_stimulus/scaling:startt,7)=Ncycle(i_dist)<2;
            
            % Execution times (Necesary for post processing of RNNs)
            exec(counter,startt:endt) = 1;
            % current cycle being performed (Necesary for post processing of RNNs)
            idx_current_cycle(counter,1:round(Ntimes/scaling)) = interp1(1:Ntimes,idx_Ncycle(idx_this_condition),linspace(1,Ntimes,round(Ntimes/scaling)));

            % identification of the conditions 
            idx_pos_trial(counter) = i_pos;
            idx_dir_trial(counter) = direction;
            idx_cycle_trial(counter) = Ncycle(i_dist);

            counter=counter+1;
        end
    end
end

%% task type is projected as two orthogonal random vectors
rng('default')
Btask=normrnd(0, 0.3,50,1);
% find an orthogonal vector to Btask
Btask2 = rotate_n_dimensional_vector(Btask,pi/2);
Btask=[Btask,Btask2]';

save(['.\Output_files\RNNs_Inputs\M1_' animal '_different.mat'],'Input','Output','idx_pos_trial','idx_cycle_trial','idx_dir_trial','exec','idx_current_cycle','Btask')
end

function Create_Inputs_same_SMA(animal,timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist)
%% Create_Inputs_same creates and saves the inputs to the RNNs that model M1 activity same control hypothesis
%% This model has 3 inputs: Start position, movement end and movement direction. Start position and movement direction are one-hot coded so that there are 5 inputs for each trial

%% 5 inputs
% 1 and 2 indicate the start of mov and start pos
% 3 indicates movement end
% 4 and 5 indicate the start of mov (same timing than 1 and 2) and
% direction

%% Parameters
Ncycle=unique(idx_dist); %[0.5, 1,2, 4,7]
NNcycle=numel(Ncycle); % 5 number of cycles
Ninputs = 5;
scaling=10; %timesteps is 10 ms
timesteps=round(5000/scaling);% 5000 ms divided by the timesteps
Ndim_output = 4; % number of dimensions of EMG data
NcondPerCycle = 4; % number of condition for each number of cycles (2 dir * 2 pos)
Ntrials = NNcycle*NcondPerCycle;
Length_stimulus = 700;% ms


%% Initialising structures
% format as in tensorflow output [trials,timesteps,Noutputs]
Output=zeros(Ntrials,timesteps,Ndim_output);

Input=zeros(Ntrials,timesteps,Ninputs);
idx_dir_trial=zeros(Ntrials,1);
idx_pos_trial=zeros(Ntrials,1);
idx_cycle_trial=zeros(Ntrials,1);
exec=zeros(Ntrials,timesteps);
idx_current_cycle=nan(Ntrials,timesteps);

counter = 1;

for i_dist = 1:NNcycle
    for direction=1:2
        for i_pos=1:2

            idx_this_condition=find(idx_dir==direction & idx_pos==i_pos & idx_dist==Ncycle(i_dist));

            % format into tensorflow output [trials,timesteps,Noutputs]
            Ntimes=numel(idx_this_condition);

            %% Define Outputs as the first 4 component of EMG
            Output(counter,1:round(Ntimes/scaling),:)=interp1(1:Ntimes,scores(idx_this_condition,1:Ndim_output),linspace(1,Ntimes,round(Ntimes/scaling)));

            % Set the final steady state
            Output(counter,round(Ntimes/scaling)+1:end,:)=repmat(Output(counter,round(Ntimes/scaling),:),1,timesteps-round(Ntimes/scaling),1);
            
            %% Define inputs
            endt = round((Ntimes-timesmov(2))/scaling); % time of movement offset
            startt = round((-timesmov(1))/scaling); % time of movement onset 

            % Inputs 1 and 2: starting position
            Input(counter, startt-round(Length_stimulus/scaling):startt,1)=i_pos-1;
            Input(counter, startt-round(Length_stimulus/scaling):startt,2)=round(-(i_pos-2)/2);

            % Input 3: temporal context input (same for all conditions)
            inputL=endt-round(100/scaling)-(startt-round(300/scaling))+1;
            
            % flat until movement onset
            Input(counter,startt-round(Length_stimulus/scaling):startt,3)=1;
            % decreasing afterwards
            Input(counter,startt-round(300/scaling):endt-round(100/scaling),3)=linspace(1,0,inputL);
            %Input(counter,endt-round(Length_stimulus/scaling)-round(100/scaling):endt-round(100/scaling),3)=1;


            % Inputs 4 and 5: movement direction
            Input(counter,startt-round(Length_stimulus/scaling):startt,4)=direction-1;
            Input(counter,startt-round(Length_stimulus/scaling):startt,5)=round(-(direction-2)/2);
            
            % Execution times (Necesary for post processing of RNNs)
            exec(counter,startt:endt) = 1;
            % current cycle being performed (Necesary for post processing of RNNs)
            idx_current_cycle(counter,1:round(Ntimes/scaling)) = interp1(1:Ntimes,idx_Ncycle(idx_this_condition),linspace(1,Ntimes,round(Ntimes/scaling)));

            % identification of the conditions 
            idx_pos_trial(counter) = i_pos;
            idx_dir_trial(counter) = direction;
            idx_cycle_trial(counter) = Ncycle(i_dist);

            counter=counter+1;
        end
    end
end

save(['.\Output_files\RNNs_Inputs\SMA_' animal '_same.mat'],'Input','Output','idx_pos_trial','idx_cycle_trial','idx_dir_trial','exec','idx_current_cycle')
end

function Create_Inputs_different_SMA(animal,timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist)
%% Create_Inputs_same creates and saves the inputs to the RNNs that model M1 activity same control hypothesis
%% This model has 3 inputs: Start position, movement end and movement direction. Start position and movement direction are one-hot coded so that there are 5 inputs for each trial

%% 5 inputs
% 1 and 2 indicate the start of mov and start pos
% 3 indicates movement end
% 4 and 5 indicate the start of mov (same timing than 1 and 2) and
% direction

%% Parameters
Ncycle=unique(idx_dist); %[0.5, 1,2, 4,7]
NNcycle=numel(Ncycle); % 5 number of cycles
Ninputs = 7;
scaling=10; %timesteps is 10 ms
timesteps=round(5000/scaling);% 5000 ms divided by the timesteps
Ndim_output = 4; % number of dimensions of EMG data
NcondPerCycle = 4; % number of condition for each number of cycles (2 dir * 2 pos)
Ntrials = NNcycle*NcondPerCycle;
Length_stimulus = 700;% ms


%% Initialising structures
% format as in tensorflow output [trials,timesteps,Noutputs]
Output=zeros(Ntrials,timesteps,Ndim_output);

Input=zeros(Ntrials,timesteps,Ninputs);
idx_dir_trial=zeros(Ntrials,1);
idx_pos_trial=zeros(Ntrials,1);
idx_cycle_trial=zeros(Ntrials,1);
exec=zeros(Ntrials,timesteps);
idx_current_cycle=nan(Ntrials,timesteps);

counter = 1;

for i_dist = 1:NNcycle
    for direction=1:2
        for i_pos=1:2

            idx_this_condition=find(idx_dir==direction & idx_pos==i_pos & idx_dist==Ncycle(i_dist));

            % format into tensorflow output [trials,timesteps,Noutputs]
            Ntimes=numel(idx_this_condition);

            %% Define Outputs as the first 4 component of EMG
            Output(counter,1:round(Ntimes/scaling),:)=interp1(1:Ntimes,scores(idx_this_condition,1:Ndim_output),linspace(1,Ntimes,round(Ntimes/scaling)));

            % Set the final steady state
            Output(counter,round(Ntimes/scaling)+1:end,:)=repmat(Output(counter,round(Ntimes/scaling),:),1,timesteps-round(Ntimes/scaling),1);
            
            %% Define inputs
            endt = round((Ntimes-timesmov(2))/scaling); % time of movement offset
            startt = round((-timesmov(1))/scaling); % time of movement onset 

            % Inputs 1 and 2: starting position
            Input(counter, startt-round(Length_stimulus/scaling):startt,1)=i_pos-1;
            Input(counter, startt-round(Length_stimulus/scaling):startt,2)=round(-(i_pos-2)/2);

            % Input 3: empty
            % Input(counter,endt-round(Length_stimulus/scaling):endt,3)=1;
            

            % Inputs 4 and 5: movement direction
            Input(counter,startt-round(Length_stimulus/scaling):startt,4)=direction-1;
            Input(counter,startt-round(Length_stimulus/scaling):startt,5)=round(-(direction-2)/2);

            % Inputs 6 and 7: indicating the type of movement (key difference with same control model)
            inputL=endt-round(100/scaling)-(startt-round(300/scaling))+1;

            % long input 
            Input(counter,startt-Length_stimulus/scaling:startt,6)=1*(Ncycle(i_dist)>=2);
            Input(counter,startt-Length_stimulus/scaling:startt,7)=1*(Ncycle(i_dist)<2);
            Input(counter,startt-round(300/scaling):endt-round(100/scaling),6)=linspace(1,0,inputL).*(Ncycle(i_dist)>=2);
            Input(counter,startt-round(300/scaling):endt-round(100/scaling),7)=linspace(1,0,inputL).*(Ncycle(i_dist)<2);

            
            % Execution times (Necesary for post processing of RNNs)
            exec(counter,startt:endt) = 1;
            % current cycle being performed (Necesary for post processing of RNNs)
            idx_current_cycle(counter,1:round(Ntimes/scaling)) = interp1(1:Ntimes,idx_Ncycle(idx_this_condition),linspace(1,Ntimes,round(Ntimes/scaling)));

            % identification of the conditions 
            idx_pos_trial(counter) = i_pos;
            idx_dir_trial(counter) = direction;
            idx_cycle_trial(counter) = Ncycle(i_dist);

            counter=counter+1;
        end
    end
end

%% task type is projected as two orthogonal random vectors
rng('default')
Btask=normrnd(0, 0.3,50,1);
% find an orthogonal vector to Btask
Btask2 = rotate_n_dimensional_vector(Btask,pi/2);
Btask=[Btask,Btask2]';

save(['.\Output_files\RNNs_Inputs\SMA_' animal '_different.mat'],'Input','Output','idx_pos_trial','idx_cycle_trial','idx_dir_trial','exec','idx_current_cycle','Btask')
end