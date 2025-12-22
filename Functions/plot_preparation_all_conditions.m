function plot_preparation_all_conditions(scores,idx_dir, idx_pos, idx_dist,column)
% PLOT_PREPARATION_CONDITIONS  Plots of preparatory neural trajectories coloured by number of cycles
%   to be performed and by directions and initial position of movement
%
%
%   PLOT_PREPARATION_CONDITIONS(scores,idx_dir, idx_pos, idx_dist,column)
%
%   INPUTS
%   ------
%   scores : [T x D] matrix
%       Neural activity matrix (RNN or neural recording), where T is the number of time samples
%       and N is the number of units or neurons.
%
%   idx_dir     : [T x 1] vector of movement or task direction labels.
%
%   idx_pos     : [T x 1] vector of position or task condition labels.
%
%   idx_dist    : [T x 1] vector indicating the total distance (Number of cycles) the animal
%       covers in the trial [0.5 1 2 4 7]
%
%   column      : scalar
%       Column index specifying where plots should be placed in the subplot
%       grid.
%
%   OUTPUT
%   ------
%   Plots of preparatory neural trajectories coloured by number of cycles
%   to be performed and by directions and initial position of movement
%
% Andrea Colins
% 22/12/2025

Ndir=2;
Npos=2;

NumberCyles=unique(idx_dist);
Ncycle=numel(NumberCyles);

colour_dir=[[115, 80, 185];[87, 204, 153]]./255;
colour_ndist=plasma(Ncycle);

%% Plot trajectories by number of cycles and then by direction 
for i_pos=1:Npos
    for i_dir=1:Ndir
        for i_cycle=1:Ncycle

            this_cond=find(idx_dist==i_cycle & idx_dir==i_dir & idx_pos==i_pos);

            % Plot trajectories by number of cycles
            subplot(2,3,3+column)
            hold on
            plot3(scores(this_cond,1),scores(this_cond,2),scores(this_cond,3),'Color',colour_ndist(i_cycle,:),'LineWidth',5/i_cycle)
            plot3(scores(this_cond(1),1),scores(this_cond(1),2),scores(this_cond(1),3),'o','Color',colour_ndist(i_cycle,:))
            
            % Plot trajectories by direction
            subplot(2,3,column)
            hold on
            plot3(scores(this_cond,1),scores(this_cond,2),scores(this_cond,3),'Color',colour_dir(i_dir,:)./sqrt(i_pos),'LineWidth',5/i_cycle)
            plot3(scores(this_cond(1),1),scores(this_cond(1),2),scores(this_cond(1),3),'o','Color',colour_dir(i_dir,:)./sqrt(i_pos))


        end
    end
end
end