function Dist_all = compare_traj_directions(scores,idx_dir,idx_pos,idx_dist)
% aim: show that trajectories of different number of cycles follow the same
% path
% compare for distances
Ndist=max(idx_dist);
Ndir=max(idx_dir);
Npos=max(idx_pos);

Nt=sum(idx_pos==1 & idx_dir==1 & idx_dist==1);
Ncomb=Ndist*(Ndist-1)./2;
total_dist=nan(Nt,Ncomb*Ndir*Npos);
counter=1;
%% across ncycles
for i_pos=1:Npos
    for i_dir=1:Ndir
        cond_ex=idx_pos==i_pos & idx_dir==i_dir;        
        total_dist(:,(Ncomb*(counter-1)+1):Ncomb*counter)=Euclidean_distance_over_time_single_cond(scores(cond_ex,:),idx_dist(cond_ex));
        counter=counter+1;
        
    end
end

%% direction
Ncomb=1;
total_dist_dir=nan(Nt,Ncomb*Ndist*Npos);
counter=1;

for i_pos=1:Npos
    for i_dist=1:Ndist
        cond_ex=idx_pos==i_pos & idx_dist==i_dist;
        total_dist_dir(:,counter)=Euclidean_distance_over_time_single_cond(scores(cond_ex,:),idx_dir(cond_ex));
        
        counter=counter+1;

    end
end

%% Position
Ncomb=1;
total_dist_pos=nan(Nt,Ncomb*Ndist*Npos);
counter=1;

for i_dir=1:Ndir
    for i_dist=1:Ndist
        cond_ex=idx_dir==i_dir & idx_dist==i_dist;
         total_dist_pos(:,counter)=Euclidean_distance_over_time_single_cond(scores(cond_ex,:),idx_pos(cond_ex));
        
        counter=counter+1;
        % plot(angles_all,'r')
        % hold on
    end
end

% maximum separation of trajectories of different initial position
max_sep_pos=max(mean(total_dist_pos,2));

% Normalise distances
Dist_cycle=mean(total_dist,2)./max_sep_pos;
Dist_dir=mean(total_dist_dir,2)./max_sep_pos;
Dist_pos=mean(total_dist_pos,2)./max_sep_pos;

Dist_all=[Dist_cycle,Dist_dir,Dist_pos];

end