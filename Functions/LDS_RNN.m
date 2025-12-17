function [epochs_all,mov_exec_all]=LDS_RNN(states,idx_dir,idx_pos,idx_dist)
timesmov=[-100 40];
t_from=timesmov(1);
ndim=3;
Ndir=2;
Npos=2;
dt=10; % delta in time between points 

threshold=50;
threshold_rec=15;

Ncycles=unique(idx_dist);
mov_exec_all=[];
Max_eigen_all=[];
idx_dir_new_all=[];
idx_pos_new_all=[];
idx_cycles_all=[];

[~,scores,~,~,explained]=pca(states);

baseline=mean(scores(1:30,:));
scores=scores-baseline;


for i_dist=1:numel(Ncycles)
    this_cond=idx_dist==Ncycles(i_dist);
    t_upto=sum(idx_dir==1 & idx_pos==1 & this_cond);
    mov_dur=(t_upto+timesmov(1)-timesmov(2))/Ncycles(i_dist); %average cycle duration


    [Max_eigen,exec,idx_dir_new,idx_pos_new,~,all_eigen,All_M]=fit_LDS_PCs(scores(this_cond,1:ndim), ...
        idx_dir(this_cond),idx_pos(this_cond),...
        explained,t_from,mov_dur*Ncycles(i_dist),threshold,threshold_rec,0,ndim^2*4);

    Max_eigen_all=[Max_eigen_all;Max_eigen];
    idx_dir_new_all=[idx_dir_new_all;idx_dir_new];
    idx_pos_new_all=[idx_pos_new_all;idx_pos_new];
    idx_cycles_all=[idx_cycles_all;ones(size(idx_pos_new,1),1)*i_dist];
    mov_exec_all=[mov_exec_all;mov_dur*ones(Ndir*Npos,1)*Ncycles(i_dist)];

end
% since dt=10 ms for these trajectories, we divide eigenvalues for dt
[~,epochs_all]=Detect_dynamics_epochs([real(Max_eigen_all),imag(Max_eigen_all)]./dt,...
    idx_dir_new_all,idx_pos_new_all,idx_cycles_all,t_from,mov_exec_all,0);

end