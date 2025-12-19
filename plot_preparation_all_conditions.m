function plot_preparation_all_conditions(states,idx_dir, idx_pos, idx_cycle,exec,column)
% PLOT_PREPARATION_CONDITIONS  Plots of preparatory neural trajectories coloured by number of cycles
%   to be performed and by directions and initial position of movement
%
%
%   PLOT_PREPARATION_CONDITIONS(states,idx_dir, idx_pos, idx_cycle,exec,column)
%
%   INPUTS
%   ------
%   states : [T x N] matrix
%       Neural activity matrix (RNN or neural recording), where T is the number of time samples
%       and N is the number of units or neurons.
%
%   idx_dir     : [T x 1] vector of movement or task direction labels.
%
%   idx_pos     : [T x 1] vector of position or task condition labels.
%
%   idx_cycle   : [T x 1] vector indicating the total distance (Number of cycles) the animal
%   covers in the trial [0.5 1 2 4 7]
%
%   exec        : [T x 1] binary vector indicating if animal is performing
%   a movement at time t (1) or not (0). 
%
%   OUTPUT
%   ------
%   Plots of preparatory neural trajectories coloured by number of cycles
%   to be performed and by directions and initial position of movement
%
% Andrea Colins Rodriguez
% 18/12/2025

Ndir=2;
Npos=2;

NumberCyles=unique(idx_cycle);
Ncycle=numel(NumberCyles);

colour_dir=[[115, 80, 185];[87, 204, 153]]./255;%[[194 165 207];[166 219 160]]./255;
colour_ndist=plasma(Ncycle);

[Ntimes,Nunits] = size(states);


counter = 1;


%% Select preparatory activity for all conditions
idxCycle_prep=nan(Ntimes,1);
idxDir_prep=nan(Ntimes,1);
idxPos_prep=nan(Ntimes,1);
FRicycle_prep_all=nan(Ntimes,Nunits);

for i_pos=1:Npos
    for i_dir=1:Ndir


        for i_cycle=1:Ncycle

            this_cond=find(idx_dir==i_dir & idx_pos==i_pos & idx_cycle==NumberCyles(i_cycle));
            this_prep=1:find(exec(this_cond)==1,1,'first');

            FRicycle=states(this_cond(this_prep),:);

            Nrows = numel(this_cond(this_prep));

            idxCycle_prep(counter:counter+Nrows-1,:)=ones(Nrows,1)*i_cycle;
            idxDir_prep(counter:counter+Nrows-1,:)=ones(Nrows,1)*i_dir;
            idxPos_prep(counter:counter+Nrows-1,:)=ones(Nrows,1)*i_pos;
            FRicycle_prep_all(counter:counter+Nrows-1,:)=FRicycle;

            counter=counter+Nrows;
        end


    end
end

idx_nan = isnan(idxCycle_prep);
idxCycle_prep(idx_nan,:) = [];
idxDir_prep(idx_nan,:) = [];
idxPos_prep(idx_nan,:) = [];
FRicycle_prep_all(idx_nan,:) = [];

%% Perform PCA only on preparation
[~,PCAPrep]=pca(FRicycle_prep_all);


%% Plot trajectories by number of cycles and then by direction 
for i_pos=1:Npos
    for i_dir=1:Ndir
        for i_cycle=1:Ncycle

            this_cond=find(idxCycle_prep==i_cycle & idxDir_prep==i_dir & idxPos_prep==i_pos);

            % Plot trajectories by number of cycles
            subplot(2,3,3+column)
            hold on
            plot3(PCAPrep(this_cond(this_prep),1),PCAPrep(this_cond(this_prep),2),PCAPrep(this_cond(this_prep),3),'Color',colour_ndist(i_cycle,:),'LineWidth',5/i_cycle)
            plot3(PCAPrep(this_cond(this_prep(1)),1),PCAPrep(this_cond(this_prep(1)),2),PCAPrep(this_cond(this_prep(1)),3),'o','Color',colour_ndist(i_cycle,:))
            
            % Plot trajectories by direction
            subplot(2,3,column)
            hold on
            plot3(PCAPrep(this_cond(this_prep),1),PCAPrep(this_cond(this_prep),2),PCAPrep(this_cond(this_prep),3),'Color',colour_dir(i_dir,:)./sqrt(i_pos),'LineWidth',5/i_cycle)
            plot3(PCAPrep(this_cond(this_prep(1)),1),PCAPrep(this_cond(this_prep(1)),2),PCAPrep(this_cond(this_prep(1)),3),'o','Color',colour_dir(i_dir,:)./sqrt(i_pos))


        end
    end
end
end