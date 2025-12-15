function [angles_all,angles_all_dir,angles_all_pos,Dist_all]=compare_traj_directions(scores,idx_dir,idx_pos,idx_dist,Npoints)
% aim: show that trajectories of different number of cycles follow the same
% path
% compare for distances
Ndist=max(idx_dist);
Ndir=max(idx_dir);
Npos=max(idx_pos);

Nt=sum(idx_pos==1 & idx_dir==1 & idx_dist==1);
Ncomb=Ndist*(Ndist-1)./2;
angles_all=nan(Nt-Npoints(2),Ncomb*Ndir*Npos);
total_dist=nan(Nt,Ncomb*Ndir*Npos);
counter=1;
%% across ncycles
for i_pos=1:Npos
    for i_dir=1:Ndir
        cond_ex=idx_pos==i_pos & idx_dir==i_dir;
        angles_all(:,(Ncomb*(counter-1)+1):Ncomb*counter)=compare_one_cond(scores(cond_ex,:),idx_dist(cond_ex),Npoints);
        
        total_dist(:,(Ncomb*(counter-1)+1):Ncomb*counter)=Euclidean_distance_over_time_single_cond(scores(cond_ex,:),idx_dist(cond_ex));
        counter=counter+1;
        
    end
end

%% direction
Ncomb=1;
angles_all_dir=nan(Nt-Npoints(2),Ncomb*Ndist*Npos);
total_dist_dir=nan(Nt,Ncomb*Ndist*Npos);
counter=1;

for i_pos=1:Npos
    for i_dist=1:Ndist
        cond_ex=idx_pos==i_pos & idx_dist==i_dist;
        angles_all_dir(:,counter)=compare_one_cond(scores(cond_ex,:),idx_dir(cond_ex),Npoints);
        total_dist_dir(:,counter)=Euclidean_distance_over_time_single_cond(scores(cond_ex,:),idx_dir(cond_ex));
        
        counter=counter+1;
        % plot(angles_all,'r')
        % hold on
    end
end

%% Position
Ncomb=1;
angles_all_pos=nan(Nt-Npoints(2),Ncomb*Ndist*Ndir);
total_dist_pos=nan(Nt,Ncomb*Ndist*Npos);
counter=1;

for i_dir=1:Ndir
    for i_dist=1:Ndist
        cond_ex=idx_dir==i_dir & idx_dist==i_dist;
        angles_all_pos(:,counter)=compare_one_cond(scores(cond_ex,:),idx_pos(cond_ex),Npoints);
        total_dist_pos(:,counter)=Euclidean_distance_over_time_single_cond(scores(cond_ex,:),idx_pos(cond_ex));
        
        counter=counter+1;
        % plot(angles_all,'r')
        % hold on
    end
end
%figure
hold on 
%plot(mean(angles_all,2),'k')
%plot(mean(angles_all_dir,2),'r')
%plot(mean(angles_all_pos,2),'g')

max_sep_pos=max(mean(total_dist_pos,2));

Dist_cycle=mean(total_dist,2)./max_sep_pos;
Dist_dir=mean(total_dist_dir,2)./max_sep_pos;
Dist_pos=mean(total_dist_pos,2)./max_sep_pos;

Dist_all=[Dist_cycle,Dist_dir,Dist_pos];
% figure
% hold on
% plot(Dist_cycle,'k')
% plot(Dist_dir,'r')
% plot(Dist_pos,'g')


end


function angles_all=compare_one_cond(scores,idx_cond,Npoints)

Ncond=max(idx_cond);
Nt=sum(idx_cond==1);
tang=nan(Nt-Npoints(2),size(scores,2),Ncond);

%compute tangent of all trajs
for i_cond=1:Ncond
    tang(:,:,i_cond)=compute_traj_dir(scores(idx_cond==i_cond,:),Npoints);

end

% compute angle between trajectories
Ncomb=Ncond*(Ncond-1)./2;
angles_all=nan(Nt-Npoints(2),Ncomb);
counter=1;
for i_cond=1:Ncond
    for j_cond=i_cond+1:Ncond
        angles_all(:,counter)=dot(tang(:,:,i_cond)',tang(:,:,j_cond)')./(vecnorm(tang(:,:,i_cond)').*vecnorm(tang(:,:,j_cond)'));
        counter=counter+1;
    end
end
% hold on
% plot(angles_all,'Color',colour_traj)
%plot(subspace(tang(:,:,1),tang(:,:,2)))

end

function tang=compute_traj_dir(scores,Npoints)
Nt=size(scores,1);
% brute force
ds=nan(Nt-Npoints(2),size(scores,2),Npoints(2));
for ipoint=Npoints(1):Npoints(2)
    ds(:,:,ipoint)=(scores(1+ipoint:end-Npoints(2)+ipoint,:)-scores(1:end-Npoints(2),:))./ipoint;
end

tang=mean(ds,3,'omitnan');

%% debugging
% scores2=scores(1:end-Npoints,:)+100*tang;
% figure
%
% for i=1:Nt-Npoints
%     plot3(scores(:,1),scores(:,2),scores(:,3))
%     hold on
%     plot3([scores(i,1) scores2(i,1)],[scores(i,2) scores2(i,2)],[scores(i,3), scores2(i,3)],'r')
%     pause(0.1)
%     cla
% end
end
