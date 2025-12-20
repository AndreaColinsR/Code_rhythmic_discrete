function [Var,Percentage_all]=dPCA_across_conditions(average_cond,idx_dir,idx_pos,idx_dist,do_plot,iplot)
% DPCA_ACROSS_CONDITIONS Performs demixed principal component analysis (dPCA) across multiple task
% conditions to quantify how population activity variance is distributed
% across task variables and their interactions.
%
% The analysis demixes neural or RNN activity into components associated
% with the following marginalisations:
%   - Movement distance (number of cycles)
%   - Movement direction
%   - Initial position
%   - Condition-independent (time-related) activity
%
%
%  DPCA_ACROSS_CONDITIONS(average_cond,idx_dir,idx_pos,idx_dist,do_plot,iplot)
%
%   INPUTS
%   ------
%   average_cond : [T x N] matrix
%       Condition-averaged neural or RNN activity, where T is the number of
%       time samples and N is the number of units or dimensions.
%
%   idx_dir      : [T x 1] vector
%       Movement direction labels.
%
%   idx_pos      : [T x 1] vector
%       Initial position condition labels.
%
%   idx_dist   : [T x 1] vector 
%       Total distance (Number of cycles) the animal
%   covers in the trial [0.5 1 2 4 7].
%
%   do_plot      : scalar (0 or 1)
%       Flag indicating whether to generate dPCA visualisations.
%
%   iplot        : scalar
%       Column index used to place plots within an existing subplot layout.
%
%   OUTPUTS
%   -------
%   Var          : vector
%       Cumulative variance explained by dPCA components.
%
%   Percentage_all : vector
%       Percentage of total variance explained by each marginalisation.
%
% Andrea Colins
% 19/12/2025

Ndir = max(idx_dir);
Npos = max(idx_pos);
Ncycle = max(idx_dist);
Nunits = size(average_cond,2);

Nt = sum(idx_dir==1 & idx_pos==1 & idx_dist==1);
FR_dPCA = nan(Nunits,Ncycle,Ndir,Npos,Nt);

for i_pos=1:Npos
    for i_dir=1:Ndir
        for i_dist=1:Ncycle
            % force to take the same number of points for all conditions
            idx=idx_dir==i_dir & idx_pos==i_pos & idx_dist==i_dist;
            % normalising the trajectories by their length
            FR_dPCA(:,i_dist,i_dir,i_pos,:) = average_cond(idx,:)';


        end
    end
end


combinedParams = {{1,[1 4]}, {2,[2 3]}, {3,[3 4]},{4}};
margNames = {'Ncycle', 'Direction','Pos','Condition-independent','D/D Interaction'};

[W,V,whichMarg] = dpca(FR_dPCA, 15, 'combinedParams', combinedParams,'lambda',1e-09);

explVar = dpca_explainedVariance(FR_dPCA, W, V, ...
    'combinedParams', combinedParams);

Var=explVar.cumulativeDPCA;
Percentage_all=round(100*(explVar.totalMarginalizedVar/explVar.totalVar));

if do_plot
    mydPCAplot(average_cond, W,whichMarg,explVar,margNames,[idx_dir,idx_dist,idx_pos],iplot)
end


end

function mydPCAplot(average_cond,W,whichMarg,explVar,margNames,conditions,column_plot)
% MYDPCAPLOT Visualises demixed principal component analysis (dPCA) results by
% projecting population activity onto demixed axes and plotting temporal
% trajectories for each task marginalisation.
%
% The function displays the first demixed component associated with each
% marginalisation and colours trajectories according to task variables.
%
% MYDPCAPLOT(average_cond,W,whichMarg,explVar,margNames,conditions,column_plot)
%
%   INPUTS
%   ------
%   average_cond : [T x N] matrix
%       Condition-averaged neural or RNN activity, where T is the number of
%       time samples and N is the number of units or neurons.
%
%   W            : [N x K] matrix
%       dPCA decoding matrix, where each column corresponds to a demixed
%       component.
%
%   whichMarg   : [K x 1] vector
%       Marginalisation index for each dPCA component.
%
%   explVar     : structure
%       Structure containing variance explained by each dPCA component and
%       marginalisation.
%
%   margNames   : cell array of strings
%       Names of marginalisations to be displayed in plot titles.
%
%   conditions  : [T x 3] matrix
%       Condition indices corresponding to direction, distance and initial
%       position.
%
%   column_plot : scalar
%       Column index specifying where plots should be placed in the subplot
%       grid.
%
%   OUTPUTS
%   -------
%   None.
%
%
% Andrea Colins
% 19/12/2025

colour_dir=[[230 165 207];[166 219 160]]./255;

%% project Neural activity onto all axes
dscore=average_cond*W;
mina=min(dscore,[],'all');
maxa=max(dscore,[],'all');
Nconditions=size(conditions,2)+1;

%Var=explVar.cumulativeDPCA;
Var_eachPC=explVar.componentVar;

idx_dir=conditions(:,1);
idx_dist=conditions(:,2);
idx_pos=conditions(:,3);


PCcond=nan(Nconditions,2);

% compute the Variance explained by the first axis of each marginalisation
for icond=1:Nconditions

    First2PC=find(whichMarg==icond,2,"first");
    if numel(First2PC)==2
        PCcond(icond,:)=First2PC;
    elseif numel(First2PC)==1
        PCcond(icond,:)=[First2PC Nconditions+1];
        Var_eachPC(Nconditions+1)=0;
        dscore(:,Nconditions+1)=0;
    else
        PCcond(icond,:)=[Nconditions+1 Nconditions+1];
        Var_eachPC(Nconditions+1)=0;
        dscore(:,Nconditions+1)=0;
    end

    subplot(Nconditions+1,3,3*(icond-1)+column_plot)
    hold on
    box off
    xlabel('Time [ms]')
    title([margNames{icond}  ', ' num2str(Var_eachPC(PCcond(icond,1)),'%.1f')])
    ylim([mina maxa])

end

Ndir=max(idx_dir);
Npos=max(idx_pos);
Ncycle=max(idx_dist);

colour_ncycles=plasma(Ncycle);

%% Plot the first component (marginal) in corresponding colour

% adapt time to timescale of the data
idx=sum(idx_dir==1 & idx_pos== 1 & idx_dist==1);
t=-idx:-1;
if idx<599
    t=-1000:10:-1;
end

for i_pos=1:Npos
    for i_dist=1:Ncycle
        for i_dir=1:Ndir
            % force to take the same number of points for all conditions
            idx=idx_dir==i_dir & idx_pos==i_pos & idx_dist==i_dist;

            for icond=1:Nconditions
                if icond==1 || icond==Nconditions
                    colour_traj=colour_ncycles(i_dist,:);

                else
                    colour_traj=colour_dir(i_dir,:)./sqrt(i_pos);
                end

                 subplot(Nconditions+1,3,3*(icond-1)+column_plot)
                plot(t,dscore(idx,PCcond(icond,1)),'Color',colour_traj,'LineWidth',2)

            end
        end
    end
end
end