function [epochs_all,mov_exec_all]=LDS_RNN(states,idx_dir,idx_pos,idx_dist)
% LDS_RNN Characterises recurrent neural network dynamics using linear
% dynamical systems.
%
%   [epochs_all,mov_exec_all] = LDS_RNN(states,idx_dir,idx_pos,idx_dist)
%
%   This function analyses recurrent neural network (RNN) activity by
%   projecting population states into a low-dimensional space and fitting
%   local linear dynamical systems (LDS) to the resulting trajectories. The
%   eigenvalues of the fitted dynamics are used to identify distinct
%   dynamical epochs across behavioural conditions.
%
%   INPUTS
%   ------
%
%   states : [T x N] matrix
%       RNN population activity, where T is the number of time bins and N is
%       the number of network units.
%
%   idx_dir : [T x 1] vector
%       Direction index for each time bin, defining movement direction
%       conditions.
%
%   idx_pos : [T x 1] vector
%       Position index for each time bin, defining movement position
%       conditions.
%
%   idx_dist: [T x 1] vector 
%       Total distance (Number of cycles) the animal
%   covers in the trial [0.5 1 2 4 7].
%
%   OUTPUTS
%   -------
%
%   epochs_all : array
%       Identified dynamical epochs derived from eigenvalue trajectories of
%       the fitted LDS models.
%
%   mov_exec_all : vector
%       Estimated movement execution duration for each condition, expressed
%       in time bins.
%
% Andrea Colins Rodriguez
% 22/12/2025

timesmov=[-100 40];
t_from=timesmov(1);
ndim=3;
Ndir=numel(unique(idx_dir));
Npos=numel(unique(idx_pos));
dt = 10; % delta in time between points 

threshold=50;

% Ncycles=unique(idx_dist);
% mov_exec_all=[];
% Max_eigen_all=[];
% idx_dir_new_all=[];
% idx_pos_new_all=[];
% idx_cycles_all=[];

[~,scores,~,~,explained]=pca(states);

baseline=mean(scores(1:30,:));
scores=scores-baseline;

Ncycles=unique(idx_dist);
mov_exec_all=nan(Ndir*Npos*numel(Ncycles),1);
Nmax = size(scores,1);
Max_eigen_all=nan(Nmax,1);
idx_dir_new_all=nan(Nmax,1);
idx_pos_new_all=nan(Nmax,1);
idx_cycles_all=nan(Nmax,1);

counter = 1;
for i_dist=1:numel(Ncycles)

    this_cond=idx_dist==Ncycles(i_dist);

    t_upto=sum(idx_dir==1 & idx_pos==1 & this_cond);

    mov_dur=(t_upto+timesmov(1)-timesmov(2))/Ncycles(i_dist); %average cycle duration


    [Max_eigen,idx_dir_new,idx_pos_new]=fit_LDS_PCs(scores(this_cond,1:ndim), ...
        idx_dir(this_cond),idx_pos(this_cond),...
        explained,threshold,ndim^2*4);

    % Max_eigen_all=[Max_eigen_all;Max_eigen];
    % idx_dir_new_all=[idx_dir_new_all;idx_dir_new];
    % idx_pos_new_all=[idx_pos_new_all;idx_pos_new];
    % idx_cycles_all=[idx_cycles_all;ones(size(idx_pos_new,1),1)*i_dist];
    % mov_exec_all=[mov_exec_all;mov_dur*ones(Ndir*Npos,1)*Ncycles(i_dist)];

    Nbins=size(idx_pos_new,1);

    Max_eigen_all(counter:counter+Nbins-1,:)=Max_eigen;
    idx_dir_new_all(counter:counter+Nbins-1,:)=idx_dir_new;
    idx_pos_new_all(counter:counter+Nbins-1,:)=idx_pos_new;
    idx_cycles_all(counter:counter+Nbins-1)=ones(Nbins,1)*i_dist;
    mov_exec_all(Ndir*Npos*(i_dist-1)+1:Ndir*Npos*i_dist,:)=mov_dur*ones(Ndir*Npos,1)*Ncycles(i_dist);
    counter = counter + Nbins+ 1;

end

idx_nan = isnan(idx_cycles_all);
idx_cycles_all(idx_nan,:) = [];
Max_eigen_all(idx_nan,:)=[];
idx_dir_new_all(idx_nan,:)=[];
idx_pos_new_all(idx_nan,:)=[];

% since dt=10 ms for these trajectories, we divide eigenvalues for dt
[~,epochs_all]=Detect_dynamics_epochs([real(Max_eigen_all),imag(Max_eigen_all)]./dt,...
    idx_dir_new_all,idx_pos_new_all,idx_cycles_all,t_from,mov_exec_all,0);

end