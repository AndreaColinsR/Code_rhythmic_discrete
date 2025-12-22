function [scorestruct,summary,Data]=jPCA_prep_exec(FR,idx_pos,idx_dir,selected_times)
% JPCA_PREP_EXEC Prepares data for jPCA analysis during movement execution.
%
% This function organises firing rate data into condition-specific structures
% suitable for jPCA, with conditions defined by initial position and movement
% direction. It is intended to test whether the planes of neural rotation
% during preparation and execution are aligned.
%
%
%   [scorestruct, summary, Data] = ...
%       JPCA_PREP_EXEC(FR, idx_pos, idx_dir, selected_times)
%
%   INPUTS
%   ------
%   FR              : [Ntimepoints x Nneurons] matrix
%       Firing rate data, where each row corresponds to a time point and
%       each column to a neural unit.
%
%   idx_pos      : [T x 1] vector
%       Initial position condition labels.
%
%   idx_dir      : [T x 1] vector
%       Movement direction labels.
%
%   selected_times  : [T x 1] or [1 x T] vector
%       Indices of time points (within each condition) to be included in the
%       jPCA analysis.
%
%   OUTPUTS
%   -------
%   scorestruct     : struct
%       Output structure returned by the jPCA routine, containing projections
%       of the data onto jPC dimensions and related scores.
%
%   summary         : struct
%       Summary information from jPCA, including jPC vectors and the fraction
%       of variance captured by each jPC. jPCs are sorted in descending order
%       of variance explained.
%
%   Data            : [Npos*Ndir x 1] struct array
%       Condition-wise data structure used as input to jPCA. Each element
%       contains:
%         - A     : [T x Nneurons] matrix of firing rates for a given
%                   position–direction condition
%         - times : [T x 1] vector of time indices corresponding
%
% Andrea Colins
% 21/12/2025


Ndim = 4;
Npos = max(idx_pos);
Ndir = max(idx_dir);

%% for jpca
jPCA_params.params=false;
jPCA_params.numPCs = Ndim;  
jPCA_params.suppressBWrosettes = true;  % these are useful sanity plots, but lets ignore them for now
jPCA_params.suppressHistograms = true;  % these are useful sanity plots, but lets ignore them for now
jPCA_params.suppressText=true;
jPCA_params.normalize=false; % Don't normalise data here because it was already normalised in a previous function.
template_struct.A = [];
template_struct.times = [];


Data = repmat(template_struct, 1,Ndir*Npos);

i_cond=1;
for i_pos=1:Npos
    for i_dir=1:Ndir
        idx=find(idx_pos==i_pos & idx_dir==i_dir);
        Data(i_cond).A=FR(idx(selected_times),:);
        Data(i_cond).times=(selected_times)';
        i_cond=i_cond+1;
    end
end

% we don't call the original version because a few functions are deprecated
[scorestruct,summary] = jPCA_new_version(Data,selected_times, jPCA_params);
[summary.varCaptEachJPC,idxtmp]=sort(summary.varCaptEachJPC,'descend');
summary.jPCs=summary.jPCs(:,idxtmp);
summary.jPCs_highD=summary.jPCs_highD(:,idxtmp);

end