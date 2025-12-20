function [idx_dir,idx_pos,idx_Ncycle,idx_dist,exec]=find_idx_conds(trials_idx,idx_conds,idx_current_cycle)
% FIND_IDX_CONDS Extracts and expands trial condition indices across time.
%
% This function maps trial-level condition labels (e.g. movement direction,
% initial position, and distance) onto time-resolved vectors aligned with
% the concatenated trial index. It also generates an execution mask that
% excludes initial and terminal time segments of each trial.
%
%
%  [idx_dir,idx_pos,idx_Ncycle,idx_dist,exec] = ...
%      FIND_IDX_CONDS(trials_idx,idx_conds,idx_current_cycle)
%
%   INPUTS
%   ------
%   trials_idx        : [Ntimepoints x 1] vector
%       Trial index for each time point, typically produced by concatenating
%       multiple trials. Each entry specifies which trial the time point
%       belongs to.
%
%   idx_conds         : [Ntrials x 3] matrix
%       Trial-level condition labels. Columns are assumed to encode:
%         - Column 1 : movement distance or condition identifier
%         - Column 2 : initial position
%         - Column 3 : movement direction
%
%   idx_current_cycle : [Ntrials x T] matrix
%       Time-resolved cycle index for each trial, specifying the current
%       movement cycle at each time point within the trial.
%
%   OUTPUTS
%   -------
%   idx_dir           : [Ntimepoints x 1] vector
%       Movement direction label assigned to each time point.
%
%   idx_pos           : [Ntimepoints x 1] vector
%       Initial position label assigned to each time point.
%
%   idx_Ncycle        : [Ntimepoints x 1] vector
%       Current movement cycle index assigned to each time point.
%
%   idx_dist          : [Ntimepoints x 1] vector
%       Movement distance (or total number of cycles) assigned to each
%       time point.
%
%   exec              : [Ntimepoints x 1] vector
%       Execution mask indicating valid execution periods. Values are set
%       to zero during the initial (first 100 samples) and terminal
%       (last 40 samples) segments of each trial, and one elsewhere.
%
% Andrea Colins Rodriguez
% 19/12/2025

Nconds=max(trials_idx);
ntimepoints=size(trials_idx,1);
idx_dir=zeros(ntimepoints,1);
idx_dist=zeros(ntimepoints,1);
idx_pos=zeros(ntimepoints,1);
idx_Ncycle=nan(ntimepoints,1);
exec=ones(ntimepoints,1);

for i_cond=1:Nconds
    this_cond=find(trials_idx==i_cond);

    idx_Ncycle(this_cond)=idx_current_cycle(i_cond,1:numel(this_cond));

    idx_dist(this_cond)=idx_conds(i_cond,1);
    idx_pos(this_cond)=idx_conds(i_cond,2);
    idx_dir(this_cond)=idx_conds(i_cond,3);

    exec(this_cond(1:100))=0;
    exec(this_cond(end-40:end))=0;

    
end

end