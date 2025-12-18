function Init_cond_t=initial_cond_helix(states,idx_dir,idx_pos,idx_dist,idx_Ncycle,exec,do_plot)


Ndist=unique(idx_dist,'stable'); % [0.5 1 2 4 7]
% Number of distances  (5)
NNdist=numel(Ndist);

%% select times of execution
exec=exec>0;
states=states(exec>0,:);

idx_dir=idx_dir(exec,:);
idx_pos=idx_pos(exec,:);
idx_dist=idx_dist(exec,:);
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

t=0:0.01:1;

counter=1;

if do_plot
    colour_dist=plasma(NNdist);
end

for i_dir=1:2
    for i_pos=1:2

        % Compute centre of each cycle within the longest rhythmic movement (7-cycles)
        for i_cycle=1:7

            this_cond=idx_dir==i_dir & idx_pos==i_pos & idx_dist==7 & idx_Ncycle==i_cycle;

            center(i_cycle,:)=mean(scores(this_cond,1:ndims));
        end

        %% fit to bezier curve
        P(1,:)=center(1,:);
        P(2,:)=center(6,:);
        P(3,:)=center(end,:);

        P1P3=[P(1,:);P(3,:)];
        P2=P(2,:);
        %% order of parameters fcn = @(coefficients,problemparameters,x,y) expression
        ftmp = @(P2) eval_bezier_min(center,P2,P1P3);
        x = fminsearch(ftmp,P2);
        P(2,:)=x;
        B=eval_bezier2(P,t);
        
        % Plot bezier curve
        if do_plot && i_dir==2 && i_pos==2
            plot3(B(:,1),B(:,2),B(:,3),'Color',[0.5 0.5 0.5])
        end


        InitC=nan(NNdist,ndims);
        idx_min=nan(NNdist,1);

        this_cond = idx_dir==i_dir & idx_pos==i_pos & idx_dist==7;

        % find the closest point 
        for i_dist=1:NNdist
            this_cond2=find(idx_dir==i_dir & idx_pos==i_pos & idx_dist==dist_sorted(i_dist));
            InitC(i_dist,:)=mean(scores(this_cond2(1:10*Nt),:),1); % use the first 100  ms of movement as the initial condition
            [~,idx_min(i_dist)]=min(pdist2(scores(this_cond,:),InitC(i_dist,:)));

            if do_plot && i_dir==2 && i_pos==2
                plot3(scores(this_cond2,1),scores(this_cond2,2),scores(this_cond2,3),'Color',colour_dist(i_dist,:))
                hold on
                plot3(scores(this_cond2(1),1),scores(this_cond2(1),2),scores(this_cond2(1),3),'.','Color',colour_dist(i_dist,:),'MarkerSize',14)

            end

        end

        % normalization by the number of timebins of the condition
        Init_cond_t(counter,:) = idx_min./sum(this_cond);

        counter=counter+1;
    end
end

end

function B=eval_bezier2(P,t_in)
Nt=numel(t_in);
B=nan(Nt,size(P,2));
for ti=1:Nt
    t=t_in(ti);
    B(ti,:)=(1-t)*(P(1,:)+t*P(2,:))+t*((1-t)*P(2,:)+t*P(3,:));
end
end

