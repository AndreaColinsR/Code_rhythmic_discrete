function Dist2Att = distance_to_attractor(states,idx_cycle,idx_Ncycle,idx_dir,idx_pos)
% DISTANCE_TO_ATTRACTOR  Compute distance of low-dimensional neural states to a limit-cycle attractor
%
%   Dist2Att = DISTANCE_TO_ATTRACTOR(states, idx_cycle, idx_Ncycle, ...
%                                    idx_dir, idx_pos)
%
%   The distance is normalised by the natural variability within the
%   attractor and the distance to the origin of the subspace
%
%   INPUTS
%   ------
%   states      : [T x N] matrix of neural activity (T time points,
%                 N neurons).
%
%   idx_cycle   : [T x 1] vector indicating the total distance (Number of cycles) the animal
%   covers in the trial [0.5 1 2 4 7]
%
%   idx_Ncycle  : [T x 1] vector indicating the current cycle number performed within a
%                 movement.
%
%   idx_dir     : [T x 1] vector of movement or task direction labels.
%
%   idx_pos     : [T x 1] vector of position or task condition labels.
%
%   OUTPUT
%   ------
%   Dist2Att    : [Ntime x Ncond x Ncycle] array containing the normalized
%                 distance to the attractor.
%                 - Ntime  : maximum number of time points (with buffer)
%                 - Ncond  : direction × position conditions
%                 - Ncycle : unique cycle types
%
% Andrea Colins
% 17/12/2025

Ndir=max(idx_dir);
Npos=max(idx_pos);

Ncycles=unique(idx_cycle);
[~,scores,~,~,explained]=pca(states);
ndims=min(find(cumsum(explained)>=80,1,'first'),6);

Ntimes=sum(idx_cycle==7  & idx_dir==1 & idx_pos==1);
Ntimes=round(Ntimes*1.1); % add 10% as buffer
Dist2Att=nan(Ntimes,4,5); % 20 conditions
counter=1;

Av_cycle=round(Ntimes/7);
Natural_dist=nan(Av_cycle,ndims,7);

for i_dir=1:Ndir
    for i_pos=1:Npos

        Rhythmic=idx_cycle==7 & idx_Ncycle>1 & idx_Ncycle<7 & idx_dir==i_dir & idx_pos==i_pos;
        %%
        % Compute the natural variability of the neural activity
        % while in the attractor

        for this_cycle=2:6
            this_idx=Rhythmic & idx_Ncycle==this_cycle;
            Scores_this_cycle=interp1(linspace(0,1,sum(this_idx)),scores(this_idx,1:ndims),linspace(0,1,Av_cycle)');
            %remove possible temporal delay using the delay between first
            %PCs
            if this_cycle>2
                if numel(Scores_this_cycle(:,1))==1
                    keyboard
                end
                [r,lag]=xcorr(Natural_dist(:,1,2),Scores_this_cycle(:,1));
                [~,minlag]=max(r);
                Scores_this_cycle=circshift(Scores_this_cycle,lag(minlag)-1,1);

            end
            Natural_dist(:,:,this_cycle)=Scores_this_cycle;
        end

        Mean_cycles=mean(Natural_dist(:,:,2:6),3);
        dist_to_mean=nan(Av_cycle,5);

        for this_cycle=2:6
            dist_to_mean(:,this_cycle-1)=min(pdist2(Natural_dist(:,:,this_cycle),Mean_cycles));
        end

        baseline_dist=mean(dist_to_mean,'all');

        for i_cycle=1:numel(Ncycles)

            Discrete=idx_cycle==Ncycles(i_cycle) & idx_dir==i_dir & idx_pos==i_pos;
            Nel=sum(Discrete);
            Mdist=pdist2(scores(Rhythmic,1:ndims),scores(Discrete,1:ndims));


            Dist2Att(1:Nel,counter,i_cycle)=min(Mdist)';
            % normalise by the distance to the origin (baseline)
            Dist2Att(:,counter,i_cycle)=max(Dist2Att(:,counter,i_cycle)-baseline_dist,0,'includenan')./(Dist2Att(1,counter,i_cycle)-baseline_dist);

        end
        counter=counter+1;
    end

end

end