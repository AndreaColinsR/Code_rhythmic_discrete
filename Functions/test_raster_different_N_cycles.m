function test_raster_different_N_cycles(animal)
%% test_raster_different_N_cycles plots the rasterplots of the three regions
% recorded. By default this functions plots the neural activity
% corresponding to the condition cycling forward starting from the bottom
% of the cycle for trials of 0.5, 2 and 7 cycles. 
%
% Input:
%
% animal: string containing the animal's name e.g. 'Drake' or 'Cousteau'
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