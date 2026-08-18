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
% 
% 17/08/2025 
% Scripts were updated to reproduce new supplementary figures

tic
%% Office
jPCA_path='C:\Users\andre\OneDrive - Universidad Adolfo Ibanez\Office computer\codes_from_papers\jPCA_ForDistribution\';
dPCA_path='C:\Users\andre\OneDrive - Universidad Adolfo Ibanez\Office computer\codes_from_papers\kobak2016';
dataset_path = 'C:\Users\andre\OneDrive - Universidad Adolfo Ibanez\Office computer\Dynamical_systems_Cortex\Data_Russo';

%% home
% tic
%  jPCA_path='C:\Users\Acer\OneDrive - Universidad Adolfo Ibanez\Office computer\codes_from_papers\jPCA_ForDistribution\';
%  dPCA_path='C:\Users\Acer\OneDrive - Universidad Adolfo Ibanez\Office computer\codes_from_papers\kobak2016';
% % 
% dataset_path = 'C:\Users\Acer\OneDrive - Universidad Adolfo Ibanez\Office computer\Dynamical_systems_Cortex\Data_Russo';
addpath('.\Functions')

%% add necessary toolboxes - not necessary for code review
addpath(genpath(jPCA_path))
addpath(genpath(dPCA_path))

% Plot supplementary figures?
plot_supp_figs.do_plot = 1; % 1 = Yes, 0 = No.

%% 0. Define PCA trajectories and save results in Output_files folder
mkdir Output_files
%create_all_output_files(dataset_path)




if plot_supp_figs.do_plot == 1
    % if plotting supplemetary figures, then create all the necessary
    % figures
    plot_supp_figs.Ouputs = {figure,figure};
    plot_supp_figs.LDS = figure;
    plot_supp_figs.Prep_M1 = figure;
    plot_supp_figs.dPCA_M1 = figure;
    plot_supp_figs.TC = figure; % temporal context
    plot_supp_figs.Prep_SMA = figure;
    plot_supp_figs.dPCA_SMA = figure;
end
%% Figure 1
animal = 'Drake';
%Plot_kinematics(animal,dataset_path)
%test_raster_different_N_cycles(animal)



%% Train neural networks here
% Create inputs and outputs for the RNNs
% Create_Inputs_RNN

%% Figure 2
region_name='M1';
output_region='EMG';
hypotheses={'same','different'};
figM1=figure;
%% Analyse neural recordings 
testing_Cortical_Data_as_RNN(region_name,figM1,plot_supp_figs)
compare_network_families(region_name,output_region,hypotheses,figM1,plot_supp_figs)


%% Video: 
RNN_name_same = 'Scores_Trained_EMG_same_Drake_14';
RNN_name_diff = 'Scores_Trained_EMG_different_Drake_12';
i_dir = 2;
i_pos = 1;
%play_video_rhythmic_discrete(region_name,animal,RNN_name_same,RNN_name_diff,i_dir,i_pos)

%% Figure 3
region_name='SMA';
output_region='M1';
hypotheses={'same','different'};
figSMA=figure;

testing_Cortical_Data_as_RNN(region_name,figSMA,plot_supp_figs)
compare_network_families(region_name,output_region,hypotheses,figSMA,plot_supp_figs)

%% Video:
animal = 'Cousteau';
RNN_name_same = 'Scores_Trained_M1_same_Cousteau_1';
RNN_name_diff = 'Scores_Trained_M1_different_Cousteau_10';
i_dir = 2;
i_pos = 1;
%play_video_rhythmic_discrete(region_name,animal,RNN_name_same,RNN_name_diff,i_dir,i_pos)


%% Supplementary figure different inputs 
if plot_supp_figs.do_plot == 1

    figure(plot_supp_figs.TC)
    Supplementary_temporal_context
    Angle_mov_type_unbiasedRNN
end
toc
end