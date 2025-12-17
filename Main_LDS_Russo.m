function Main_LDS_Russo(animal,do_plot)

%% Calculate the eigenvalues of the piecewise LDS during execution (in the segment defined by Ames)
timesmov=[-1000 400]; %Ames used 1250:3250 => cycles 2 to 5
threshold=50;
threshold_rec=15;
ms=1000;
Ndir=2;
Npos=2;
t_from=timesmov(1)/ms;

Max_eigen_all=[];
idx_dir_new_all=[];
idx_pos_new_all=[];
idx_cycles_all=[];
mov_exec_all=[];

load(['.\Output_files\scores_' animal '_M1.mat'],'scores','idx_dir','idx_pos','idx_dist','explained')

ndim=find(cumsum(explained)>=threshold,1,'first');
Ncycles=unique(idx_dist);

for i_dist=1:numel(Ncycles)

    t_upto=sum(idx_dir==1 & idx_pos==1 & idx_dist==Ncycles(i_dist))/ms;
    mov_dur=(t_upto+timesmov(1)/ms-timesmov(2)/ms)*ms/Ncycles(i_dist); %average cycle duration


    [Max_eigen,~,idx_dir_new,idx_pos_new]=fit_LDS_PCs(scores(idx_dist==Ncycles(i_dist),:), ...
        idx_dir(idx_dist==Ncycles(i_dist)),idx_pos(idx_dist==Ncycles(i_dist)),explained,t_from,mov_dur*Ncycles(i_dist)/ms,threshold,threshold_rec,0,ndim^2*20);

    Max_eigen_all=[Max_eigen_all;Max_eigen];
    idx_dir_new_all=[idx_dir_new_all;idx_dir_new];
    idx_pos_new_all=[idx_pos_new_all;idx_pos_new];
    idx_cycles_all=[idx_cycles_all;ones(size(idx_pos_new,1),1)*i_dist];
    mov_exec_all=[mov_exec_all;mov_dur*ones(Ndir*Npos,1)*Ncycles(i_dist)];

end

Detect_dynamics_epochs([real(Max_eigen_all),imag(Max_eigen_all)],idx_dir_new_all,idx_pos_new_all,idx_cycles_all,timesmov(1),mov_exec_all,do_plot);

end