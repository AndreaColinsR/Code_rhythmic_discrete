function Init_cond_t = initial_cond_helix(states,idx_dir,idx_pos,idx_cycle,idx_Ncycle,exec,do_plot)
% INITIAL_COND_HELIX  Computes the position of a low dimensional neural
% state along the helix in SMA at movement onset 
%
%   Init_cond_t = INITIAL_COND_HELIX(states, idx_dir, idx_pos, ...
%                                    idx_dist, idx_Ncycle,exec,do_plot)
%
%
%   INPUTS
%   ------
%   states      : [T x N] matrix of neural activity (T time points,
%                 N neurons).
%
%   idx_dir     : [T x 1] vector of movement or task direction labels.
%
%   idx_pos     : [T x 1] vector of position or task condition labels.
%
%   idx_cycle   : [T x 1] vector indicating the total distance (Number of cycles) the animal
%   covers in the trial [0.5 1 2 4 7]
%
%   idx_Ncycle  : [T x 1] vector indicating the current cycle number performed within a
%                 movement.
%
%   do_plot : logical scalar
%       If do_plot equals 1, a 3-D PCA visualization of discrete and rhythmic trajectories
%       is generated.
%   OUTPUT
%   ------
%   Init_cond_t   : [Ncond (4) x Ncycle (5)] array containing the normalized
%                 distance to the attractor.
%
% Andrea Colins
% 18/12/2025

Ndist=unique(idx_cycle,'stable'); % [0.5 1 2 4 7]
% Number of distances  (5)
NNdist=numel(Ndist);

%% select times of execution
exec=exec>0;
states=states(exec>0,:);

idx_dir=idx_dir(exec,:);
idx_pos=idx_pos(exec,:);
idx_cycle=idx_cycle(exec,:);
idx_Ncycle=idx_Ncycle(exec,:);

Init_cond_t=nan(2*2,NNdist);

%% Execution subspace
[~,scores]=pca(states);

Nt=1;
if size(states,1)>10000
    % take more points if the trajectory correspond to dt = 1 ms (cortical trajectories)
    Nt=10;
end

ndims=6;
scores=scores(:,1:ndims);
center=nan(7,ndims);

[dist_sorted,~]=sort(Ndist,'ascend');

s=0:0.01:1;
Ns=numel(s);

counter=1;

if do_plot
    colour_dist=plasma(NNdist);
end

for i_dir=1:2
    for i_pos=1:2

        % Compute centre of each cycle within the longest rhythmic movement (7-cycles)
        for i_cycle=1:7

            this_cond=idx_dir==i_dir & idx_pos==i_pos & idx_cycle==7 & idx_Ncycle==i_cycle;

            center(i_cycle,:)=mean(scores(this_cond,1:ndims));
        end

        %% fit to bezier curve
        P(1,:)=center(1,:);
        P(2,:)=center(6,:); % Initial guess of the middle point before the fitting.
        P(3,:)=center(end,:);

        P1P3=[P(1,:);P(3,:)];
        P2=P(2,:);

        %% order of parameters fcn = @(coefficients,problemparameters,x,y) expression
        ftmp = @(P2) eval_bezier_min(center,P2,P1P3);
        middle_point = fminsearch(ftmp,P2);
        P(2,:) = middle_point;

        % evaluate
        B = eval_bezier(P,s);
        
        % Plot bezier curve
        if do_plot && i_dir==2 && i_pos==2
            plot3(B(:,1),B(:,2),B(:,3),'Color',[0.5 0.5 0.5])
        end


        InitC=nan(NNdist,ndims);
        idx_min=nan(NNdist,1);

        %% find the closest point to the Bezier curve

        for i_dist=1:NNdist

            this_cond2=find(idx_dir==i_dir & idx_pos==i_pos & idx_cycle==dist_sorted(i_dist));

            % compute the average position of the initial condition 
            InitC(i_dist,:)=mean(scores(this_cond2(1:5*Nt),:),1); % use the first 100  ms of movement as the initial condition
            [~,idx_min(i_dist)]=min(pdist2(B,InitC(i_dist,:))); % idx_min is s at which the in the Bezier curve is closest to the point

            if do_plot && i_dir==2 && i_pos==2
                plot3(scores(this_cond2,1),scores(this_cond2,2),scores(this_cond2,3),'Color',colour_dist(i_dist,:))
                hold on
                plot3(scores(this_cond2(1),1),scores(this_cond2(1),2),scores(this_cond2(1),3),'.','Color',colour_dist(i_dist,:),'MarkerSize',14)

            end

        end

        Init_cond_t(counter,:) = idx_min/Ns;
        
        counter=counter+1;
    end
end

end

