function Main_reproduce_figures
%% Main_reproduce_figures reproduces the first figure of the manuscript
% for this, download a copy of the data from the Mendeley repository.
%
% This code will create a new folder named "Output_files" to save
% intermediate results
%
%
% 04/09/2025
% Andrea Colins Rodriguez

%% Office
%jPCA_path='C:\Users\andrea.colins\OneDrive - Universidad Adolfo Ibanez\Office computer\codes_from_papers\jPCA_ForDistribution\';
%dPCA_path='C:\Users\andrea.colins\OneDrive - Universidad Adolfo Ibanez\Office computer\codes_from_papers\kobak2016';

%% home

jPCA_path='C:\Users\Acer\OneDrive - Universidad Adolfo Ibanez\Office computer\codes_from_papers\jPCA_ForDistribution\';
dPCA_path='C:\Users\Acer\OneDrive - Universidad Adolfo Ibanez\Office computer\codes_from_papers\kobak2016';

dataset_path = 'C:\Users\Acer\OneDrive - Universidad Adolfo Ibanez\Office computer\Dynamical_systems_Cortex\Data_Russo';

plot_supp = 0; % 1 to plot supplementary figures, 0 otherwise


%% add necessary toolboxes - not necessary for code review
%%addpath(genpath('C:\Users\andrea.colins\OneDrive - Universidad Adolfo Ibanez\Office computer\codes_from_papers\jPCA_ForDistribution'))
%%addpath(genpath('C:\Users\andrea.colins\OneDrive - Universidad Adolfo Ibanez\Office computer\codes_from_papers\kobak2016'))



%% 0. Define PCA trajectories and save results in Output_files folder
mkdir Output_files
%create_all_output_files(dataset_path)

%% Figure 1
animal = 'Drake';
%Plot_kinematics(animal,dataset_path)
%test_raster_different_N_cycles(animal)



%% Train neural networks here
% Create inputs and outputs for the RNNs


%% Figure 2
region_name='M1';
%compare_network_families(region_name,plot_supp)


%% Video: 
RNN_name_same = 'Scores_Trained_EMG_Hyp_continuousDrake_14_SLen300v2';
RNN_name_diff = 'Scores_Trained_EMG_Hyp_separateDrake_12_Orth_start';
i_dir = 2;
i_pos = 1;
%play_video_rhythmic_discrete(region_name,animal,RNN_name_same,RNN_name_diff,i_dir,i_pos)

%% Figure 3
region_name='SMA';
% compare_network_families(region_name,plot_supp)

%% Video:
animal = 'Cousteau';
RNN_name_same = 'Scores_Trained_M1_Hyp_continuousCousteau_1';
RNN_name_diff = 'Scores_Trained_M1_Hyp_separateCousteau_10_Orth_start';
i_dir = 2;
i_pos = 1;
%play_video_rhythmic_discrete(region_name,animal,RNN_name_same,RNN_name_diff,i_dir,i_pos)
end