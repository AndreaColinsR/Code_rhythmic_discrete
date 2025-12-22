function test_raster_different_N_cycles(animal)
% TEST_RASTER_DIFFERENT_N_CYCLES Plots raster plots for different cycle counts for one recording.
%
% This function generates raster (heatmap) plots of neural activity recorded
% from three regions (SMA, M1, and EMG) for trials with different numbers of
% movement cycles. By default, it visualises activity for forward movements
% starting from the top position, comparing selected cycle counts.Neural 
% units are sorted according to the timing of their peak activity in
% the highest cycle condition, allowing consistent ordering across plots.
%
%
%   TEST_RASTER_DIFFERENT_N_CYCLES(animal)
%
%   INPUTS
%   ------
%   animal : char or string
%       Name of the animal used to identify the corresponding data files
%       (e.g. 'Drake' or 'Cousteau').
%
%   OUTPUTS
%   -------
%   None
%       The function produces a figure containing multiple raster plots
%       arranged by brain region and number of cycles.
%
%   DETAILS
%   -------
%   - Data are loaded from precomputed score files corresponding to each
%     region (SMA, M1, EMG).
%   - Raster plots are generated for a subset of available cycle counts
%     (e.g. 0.5, 2, and 7 cycles), plotted in descending order.
%   - Neurons are sorted based on the time of peak firing rate within a
%     specified window for the largest cycle condition, and this ordering
%     is reused across conditions.
%   - Time is displayed relative to movement onset.
%
% 13/01/2025
% Andrea Colins Rodriguez

counter=10;
i_dir=2; % forward 
i_pos=1; % top
region_name={'SMA','M1','EMG'};

for i_region = 1:size(region_name,2)
load(['.\Output_files\scores_' animal '_'  region_name{i_region} '.mat'],'FR','idx_dir','idx_pos','idx_dist')

Ndist=unique(idx_dist); % 0.5,1,2,4,7
NNdist=numel(Ndist);

Nunits=size(FR,2);

for i_dist=NNdist:-2:1

    cond=idx_dist==Ndist(i_dist) & idx_dir==i_dir & idx_pos==i_pos;

    FRsmall=FR(cond,:);

    % sort rows of each heatmap by timing of FR of the neurons in the
    % condition Ndist = 7. 
    if i_dist==NNdist && i_pos==1 && i_dir==2

        [~,idxmax]=max(FRsmall(1000:1500,:));
        [~,idxsort]=sort(idxmax);

    end

    subplot(6,3,counter)
    
    imagesc([-1000 size(FRsmall,1)-1000],[1 Nunits],FRsmall(:,idxsort)')
    box off
    title([ 'N cycles = ' num2str(Ndist(i_dist))])
    ylabel(region_name{i_region})
    
    counter=counter+1;
    if i_region==3
    xlabel('Time to movement onset [ms]')
    end
end
end

colormap(flipud(colormap(gray)))
end