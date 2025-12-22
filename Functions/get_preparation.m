function prepdata = get_preparation(states,idx_dir, idx_pos, idx_cycle,exec)
% GET_PREPARATION  Extract preparatory neural trajectories for all
% conditions
%
%
%   GET_PREPARATION(states,idx_dir, idx_pos, idx_cycle,exec)
%
%   INPUTS
%   ------
%   states      : [T x N] matrix
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
            this_prep=1:(find(exec(this_cond)==1,1,'first')-1);

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
[coeff_prep,score_prep]=pca(FRicycle_prep_all);


prepdata.FR=FRicycle_prep_all;
prepdata.ndist=idxCycle_prep;
prepdata.ndir=idxDir_prep;
prepdata.npos=idxPos_prep;
prepdata.scores=score_prep;
prepdata.coeff=coeff_prep;

end




