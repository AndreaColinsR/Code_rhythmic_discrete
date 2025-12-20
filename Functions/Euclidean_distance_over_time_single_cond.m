function [total_dist,comb]=Euclidean_distance_over_time_single_cond(score,idx_condition)
% EUCLIDEAN_DISTANCE_OVER_TIME_SINGLE_COND Computes pairwise Euclidean distances over time between population
% trajectories corresponding to different levels of a single task
% condition. For each pair of condition values, the function calculates the Euclidean
% distance between low-dimensional activity trajectories, ensuring that
% all trajectories are compared over an equal number of time samples.
%
% EUCLIDEAN_DISTANCE_OVER_TIME_SINGLE_COND(score,idx_condition)
% 
%
%   INPUTS
%   ------
%   score       : [T x D] matrix
%       Low-dimensional neural or RNN trajectory (e.g. PCA scores), where T
%       is the number of time samples and D is the number of dimensions.
%
%   idx_condition : [T x 1] vector
%       Vector of condition labels defining the task variable to be
%       compared (e.g. direction, position or distance).
%
%   OUTPUTS
%   -------
%   total_dist  : [T_min x C] matrix
%       Euclidean distance over time between all pairs of condition
%       trajectories, where T_min is the minimum number of samples shared
%       across conditions and C is the number of unique condition pairs.
%
%   comb        : [C x 2] matrix
%       Index of condition pairs used to compute distances. Each row
%       specifies the two condition labels compared.
%
%
% Andrea Colins 
% 19/12/2025


% If trajectories have different lenghts, then cut all of them short to match the shortest one
Ncond=max(idx_condition);
Nsamples=nan(Ncond,1);

for i_cond=1:Ncond

    Nsamples(i_cond)=sum(idx_condition==i_cond);
end

Nsamples=min(Nsamples);

Ncomb=sum(1:Ncond-1);
total_dist=nan(Nsamples,Ncomb);
comb=nan(Ncomb,2);
counter=1;

% Compute distances for between all pairs of trajectories
for i_cond=1:Ncond

    idx1=find(idx_condition==i_cond);
 
    for j_cond=i_cond+1:Ncond
        idx2=find(idx_condition==j_cond);
        
        total_dist(:,counter)=sqrt(sum((score(idx1(1:Nsamples),:)-score(idx2(1:Nsamples),:)).^2,2));
        comb(counter,:)=[i_cond,j_cond];
        counter=counter+1;
      
    end
end
end