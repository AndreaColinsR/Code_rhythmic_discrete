function create_all_output_files(dataset_path)
%% CREATE_ALL_FILES generates a file for each animal and cortical area recorded.
% 
%   CREATE_ALL_FILES(dataset_path)
%
%   INPUTS
%   ------
%   dataset_path     : Character array
%                     path where the dataset is stored
%
%   OUTPUT
%   ------
%   Files for each animal and physiological recording (SMA, M1 and EMG).
%   Each file contains:
%
%%  Behavioural information
%
%       exec: array of N_Times elementes denoting if in a given timebin the animal was executing a movement (1) or not (0)
%       idx_dir: array of N_Times elementes denoting if in a given timebin corresponded to a trial when the animal cycled forward (2) or backward (1)
%       idx_dist: array of N_Times elementes indicating the pedalling distance of the trial (0.5, 1, 2, 4, 7)
%       idx_Ncycle: array of N_Times elementes indicating the number of the current cycle performed (1 to 7, Nan if not executing)
%       idx_pos: array of N_Times elementes indicating the starting position of the trial: bottom start (1), top start (2)
%
%% Neural data
%
%       FR: Matrix [N_times x N_neurons] Average firing rate of each neuron for
%       each trial type. FR published was already filtered with a gaussian filter
%       (sigma = 25).
%
%       scores: Matrix [N_times x 10] of neural activity projected onto the first 10 Principal Components (PCs)
%
%       explained: Array [N_neurons] indicating the variance explained by each
%       PC
%
% Both behavioural and neural data are selected from 1000 ms before
% movement onset to 400 ms after movement offset
% 
% Andrea Colins Rodriguez
% 19/12/2025 

animal={'Cousteau','Drake'};
region_name={'SMA','M1','EMG'};

timesmov = [-1000 400]; % select times from 1s before movement onset up to 400 ms after movement offset
ndims = 10; %PCs dimensions to save

for i_animal=1:numel(animal)

    for i_region=1:numel(region_name)

        [scores,explained,idx_pos,idx_dir,idx_dist,baseline,FR,idx_Ncycle,exec]=extract_trajectories_all(animal{i_animal},region_name{i_region},timesmov,dataset_path);
        % impose that the origin of the subspace correspond to the average neural
        % activity of the population 1.5 s before the movement start
        scores=scores-baseline;
        scores=scores(:,1:ndims);

        %% save into a file
        save(['.\Output_files\scores_' animal{i_animal} '_' region_name{i_region} '.mat'],'scores','explained', 'idx_dir','idx_Ncycle','idx_pos','idx_dist','FR','exec')

    end
end
end