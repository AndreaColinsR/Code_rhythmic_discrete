function Create_Inputs_RNN
% CREATE_INPUTS_RNN  Generate RNN input datasets to train RNN models using
% Tensorflow
%
%   CREATE_INPUTS_RNN generates and saves all input and output structures
%   required to train recurrent neural networks (RNNs) modelling motor
%   cortical activity. The function prepares datasets for both primary
%   motor cortex (M1) and supplementary motor area (SMA) models, under
%   different control hypotheses.
%
%
% INPUT DATA FILES
% ----------------
% The function expects the following MAT-files to exist:
%
%   .\Output_files\scores_<animal>_<output_region>.mat
%
% This files is created by the function CREATE_ALL_OUTPUT_FILES.
% Each file must contain:
%   - scores (PCA neural trajectories)
%   - idx_dir
%   - idx_Ncycle
%   - idx_pos
%   - idx_dist
%
% OUTPUT
% -------------------------------------------------------------------------
% This function does not return variables to the MATLAB workspace.
% All outputs are saved to disk by the corresponding input-generation
% functions in:
%
%   .\Output_files\RNNs_Inputs\
%
%
%
% Andrea Colins Rodriguez
% 19/12/2025

animal={'Cousteau','Drake'}; 
timesmov=[-1000 400];

%% For M1 models
output_region = 'EMG';

for i_animal=1:numel(animal)
    load(['.\Output_files\scores_' animal{i_animal} '_' output_region '.mat'],'scores', 'idx_dir','idx_Ncycle','idx_pos','idx_dist')

    Create_Inputs_same_M1(animal{i_animal},timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist);
    Create_Inputs_different_M1(animal{i_animal},timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist);

end

%% For SMA models
output_region = 'M1';

for i_animal=1:numel(animal)
    load(['.\Output_files\scores_' animal{i_animal} '_' output_region '.mat'],'scores', 'idx_dir','idx_Ncycle','idx_pos','idx_dist')

    Create_Inputs_same_SMA(animal{i_animal},timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist);
    Create_Inputs_different_SMA(animal{i_animal},timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist);

end

end

function Create_Inputs_same_M1(animal,timesmov,scores,idx_dir,idx_pos,idx_Ncycle,idx_dist)
% CREATE_INPUTS_SAME_M1  Create and save inputs for RNNs modelling M1
% activity under the "same control" hypothesis
%
%   The model uses five input channels describing movement-related events:
%     1–2 : Starting position (one-hot encoded)
%     3   : Movement end
%     4–5 : Movement direction (one-hot encoded)
%
%   Start position and movement direction are encoded using two channels
%   each, with activity occurring around movement onset. Movement end is
%   encoded using a single channel with activity around movement offset.
%
%   The function interpolates EMG data to a fixed number of time steps,
%   constructs corresponding input signals, and saves the resulting
%   variables to a MAT-file.
%
%   INPUTS
%   ------
%   animal          : Character array
%                     Identifier of the animal (used for file naming)
%
%   timesmov        : Numeric vector [1 x 2]
%                     Movement onset and offset times (in ms)
%
%   scores          : [T x D] matrix
%                     EMG neural trajectories matrix where T is the number of time samples
%                     and D is the number of dimensions; the first four dimensions are used
%                     as output dimensions
%
%   idx_dir         : [T x 1] vector of movement or task direction labels.
%
%   idx_pos         : [T x 1] vector of position or task condition labels.
%
%   idx_Ncycle      : [T x 1] vector indicating the current cycle number performed within a
%                       movement.
%
%   idx_dist        : [T x 1] vector indicating the total distance (Number of cycles) the animal
%                    covers in the trial [0.5 1 2 4 7]
%
%
%   OUTPUT
%   ------
%   Input           : [trials (20) × timesteps × 5 (Channels input)] Numeric array
%                     RNN input tensor of size
%                     
%
%   Output          : [trials × timesteps × 4 (Dimensions)] Numeric array
%                     RNN output tensor (EMG) of size
%                     
%
%   idx_pos_trial   : Numeric vector
%                     Starting position associated with each trial
%
%   idx_dir_trial   : Numeric vector
%                     Movement direction associated with each trial
%
%   idx_cycle_trial : Numeric vector
%                     Cycle condition associated with each trial
%
%   exec            : Numeric array
%                     Binary matrix indicating movement execution period
%                     for each trial
%
%   idx_current_cycle : Numeric array
%                       Interpolated cycle index over time for each trial
%
%   The outputs are saved to:
%     .\Output_files\RNNs_Inputs\M1_<animal>_same.mat
%
%
% Andrea Colins Rodriguez
% 19/12/2025

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
% CREATE_INPUTS_DIFFERENT_M1  Create and save inputs for RNNs modelling M1
% activity under the "different control" hypothesis
%
%   The model uses five input channels describing movement-related events:
%     1–2 : Starting position (one-hot encoded)
%     3   : Movement end
%     4–5 : Movement direction (one-hot encoded)
%     6–7 : Movement type (rhythmic or discrete, one-hot encoded)
%
%   Start position, movement direction and movement type are encoded using two channels
%   each, with activity occurring around movement onset. Movement end is
%   encoded using a single channel with activity around movement offset.
%
%   The function interpolates EMG data to a fixed number of time steps,
%   constructs corresponding input signals, and saves the resulting
%   variables to a MAT-file.
%
%   INPUTS
%   ------
%   animal          : Character array
%                     Identifier of the animal (used for file naming)
%
%   timesmov        : Numeric vector [1 x 2]
%                     Movement onset and offset times (in ms)
%
%   scores          : [T x D] matrix
%                     EMG neural trajectories matrix where T is the number of time samples
%                     and D is the number of dimensions; the first four dimensions are used
%                     as output dimensions
%
%   idx_dir         : [T x 1] vector of movement or task direction labels.
%
%   idx_pos         : [T x 1] vector of position or task condition labels.
%
%   idx_Ncycle      : [T x 1] vector indicating the current cycle number performed within a
%                       movement.
%
%   idx_dist        : [T x 1] vector indicating the total distance (Number of cycles) the animal
%                    covers in the trial [0.5 1 2 4 7]
%
%
%   OUTPUT
%   ------
%   Input           : [trials (20) × timesteps × 7 (Channels input)] Numeric array
%                     RNN input tensor of size
%                     
%
%   Output          : [trials × timesteps × 4 (Dimensions)] Numeric array
%                     RNN output tensor (EMG) of size
%                     
%
%   idx_pos_trial   : Numeric vector
%                     Starting position associated with each trial
%
%   idx_dir_trial   : Numeric vector
%                     Movement direction associated with each trial
%
%   idx_cycle_trial : Numeric vector
%                     Cycle condition associated with each trial
%
%   exec            : Numeric array
%                     Binary matrix indicating movement execution period
%                     for each trial
%
%   idx_current_cycle : Numeric array
%                       Interpolated cycle index over time for each trial
% 
%   Btask           : [2 x Neurons] Numeric array
%                     Two orthogonal vectors that can be used as fixed
%                     weights for the movement type inputs. 
%
%   The outputs are saved to:
%     .\Output_files\RNNs_Inputs\M1_<animal>_different.mat
%
%
% Andrea Colins Rodriguez
% 19/12/2025

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
% CREATE_INPUTS_SAME_SMA  Create and save inputs for RNNs modelling SMA
% activity under the "same control" hypothesis
%
%   The model uses five input channels describing movement-related events:
%     1–2 : Starting position (one-hot encoded)
%     3   : Temporal context
%     4–5 : Movement direction (one-hot encoded)
%
%   Start position and movement direction are encoded using two channels
%   each, with activity occurring around movement onset. Temporal context is
%   encoded using a single channel with activity that is stable during preparation and decreasing during execution.
%
%   The function interpolates M1 data to a fixed number of time steps,
%   constructs corresponding input signals, and saves the resulting
%   variables to a MAT-file.
%
%   INPUTS
%   ------
%   animal          : Character array
%                     Identifier of the animal (used for file naming)
%
%   timesmov        : Numeric vector [1 x 2]
%                     Movement onset and offset times (in ms)
%
%   scores          : [T x D] matrix
%                     M1 neural trajectories matrix where T is the number of time samples
%                     and D is the number of dimensions; the first four dimensions are used
%                     as output dimensions
%
%   idx_dir         : [T x 1] vector of movement or task direction labels.
%
%   idx_pos         : [T x 1] vector of position or task condition labels.
%
%   idx_Ncycle      : [T x 1] vector indicating the current cycle number performed within a
%                       movement.
%
%   idx_dist        : [T x 1] vector indicating the total distance (Number of cycles) the animal
%                    covers in the trial [0.5 1 2 4 7]
%
%
%   OUTPUT
%   ------
%   Input           : [trials (20) × timesteps × 5 (Channels input)] Numeric array
%                     RNN input tensor of size
%                     
%
%   Output          : [trials × timesteps × 4 (Dimensions)] Numeric array
%                     RNN output tensor (M1) of size
%                     
%
%   idx_pos_trial   : Numeric vector
%                     Starting position associated with each trial
%
%   idx_dir_trial   : Numeric vector
%                     Movement direction associated with each trial
%
%   idx_cycle_trial : Numeric vector
%                     Cycle condition associated with each trial
%
%   exec            : Numeric array
%                     Binary matrix indicating movement execution period
%                     for each trial
%
%   idx_current_cycle : Numeric array
%                       Interpolated cycle index over time for each trial
%
%   The outputs are saved to:
%     .\Output_files\RNNs_Inputs\SMA_<animal>_same.mat
%
%
% Andrea Colins Rodriguez
% 19/12/2025

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
% CREATE_INPUTS_DIFFERENT_SMA  Create and save inputs for RNNs modelling SMA
% activity under the "same control" hypothesis
%
%   The model uses five input channels describing movement-related events:
%     1–2 : Starting position (one-hot encoded)
%     3   : Empty (to match M1 model)
%     4–5 : Movement direction (one-hot encoded)
%     6–7 : Movement type, temporal context (rhythmic or discrete, one-hot encoded)
%
%   Start position, movement direction are encoded using two channels
%   each, with activity occurring around movement onset. Movement type is
%   encoded using a two channels with stable activity during preparation and 
%   decreasing activity during execution. This way, these channels also indicate temporal context.
%
%   The function interpolates M1 data to a fixed number of time steps,
%   constructs corresponding input signals, and saves the resulting
%   variables to a MAT-file.
%
%   INPUTS
%   ------
%   animal          : Character array
%                     Identifier of the animal (used for file naming)
%
%   timesmov        : Numeric vector [1 x 2]
%                     Movement onset and offset times (in ms)
%
%   scores          : [T x D] matrix
%                     M1 neural trajectories matrix where T is the number of time samples
%                     and D is the number of dimensions; the first four dimensions are used
%                     as output dimensions
%
%   idx_dir         : [T x 1] vector of movement or task direction labels.
%
%   idx_pos         : [T x 1] vector of position or task condition labels.
%
%   idx_Ncycle      : [T x 1] vector indicating the current cycle number performed within a
%                       movement.
%
%   idx_dist        : [T x 1] vector indicating the total distance (Number of cycles) the animal
%                    covers in the trial [0.5 1 2 4 7]
%
%
%   OUTPUT
%   ------
%   Input           : [trials (20) × timesteps × 5 (Channels input)] Numeric array
%                     RNN input tensor of size
%                     
%
%   Output          : [trials × timesteps × 4 (Dimensions)] Numeric array
%                     RNN output tensor (M1) of size
%                     
%
%   idx_pos_trial   : Numeric vector
%                     Starting position associated with each trial
%
%   idx_dir_trial   : Numeric vector
%                     Movement direction associated with each trial
%
%   idx_cycle_trial : Numeric vector
%                     Cycle condition associated with each trial
%
%   exec            : Numeric array
%                     Binary matrix indicating movement execution period
%                     for each trial
%
%   idx_current_cycle : Numeric array
%                       Interpolated cycle index over time for each trial
%
%   Btask           : [2 x Neurons] Numeric array
%                     Two orthogonal vectors that can be used as fixed
%                     weights for the movement type inputs. 
%
%   The outputs are saved to:
%     .\Output_files\RNNs_Inputs\SMA_<animal>_different.mat
%
%
% Andrea Colins Rodriguez
% 19/12/2025

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