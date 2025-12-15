function [Var,Percentage_all]=dPCA_across_conditions(average_cond,idx_dir,idx_pos,idx_dist,do_plot,iplot)
% aim: perform dPCA across conditions

%% dpca bit
%- 1 duration
%- 2 direction
%- 3 initial position
%- 4 time
% normalise FR

Ndir=max(idx_dir);
Npos=max(idx_pos);
Ncycle=max(idx_dist);

for i_pos=1:Npos
    for i_dir=1:Ndir
        for i_dist=1:Ncycle
            % force to take the same number of points for all conditions
            idx=idx_dir==i_dir & idx_pos==i_pos & idx_dist==i_dist;
            % normalising the trajectories by their length
            FR_dPCA(:,i_dist,i_dir,i_pos,:)=average_cond(idx,:)';


        end
    end
end

%combinedParams = {{1,[1 3]}, {2,[2 4]}, {3,[3 4]},{4}};
combinedParams = {{1,[1 4]}, {2,[2 3]}, {3,[3 4]},{4}};
margNames = {'Ncycle', 'Direction','Pos','Condition-independent','D/D Interaction'};
%margColours =  [0.5 0.5 0.5; 1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1; 0 0 0];

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