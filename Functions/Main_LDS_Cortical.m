function Main_LDS_Cortical(animal,timesmov,do_plot)
% MAIN_LDS_CORTICAL Estimates piecewise linear dynamical systems for
% cortical population activity across task conditions.
%
%   Main_LDS_Cortical(animal,timesmov,do_plot)
%
%   The function loads low-dimensional cortical activity for a given
%   animal, fits local linear dynamical systems (LDS) to neural trajectories
%   across movement conditions, and extracts eigenvalue-based dynamical
%   features. These features are subsequently used to identify and
%   optionally visualise distinct dynamical epochs during movement.
%
%   INPUTS
%   ------
%
%   animal : character vector or string
%       Identifier for the animal/session used to load precomputed neural
%       scores from disk.
%
%   timesmov : [1 x 2] vector
%       Time markers (in milliseconds) defining movement-related epochs,
%       typically corresponding to movement onset and offset.
%
%   do_plot : logical
%       Flag indicating whether detected dynamical epochs should be plotted.
%
%   OUTPUTS
%   -------
%
%   None
%
%   NOTES
%   -----
%   The function loads low-dimensional M1 activity (e.g. PCA scores) and
%   associated task indices, determines the number of dimensions required
%   to explain a specified proportion of variance, and applies
%   fit_LDS_PCs to each movement cycle condition. 
%
% Andrea Colins Rodriguez
% 22/12/2025


%% Calculate the eigenvalues of the piecewise LDS 
threshold=50;
ms=1000;
Ndir=2;
Npos=2;

load(['.\Output_files\scores_' animal '_M1.mat'],'scores','idx_dir','idx_pos','idx_dist','explained')


ndim=find(cumsum(explained)>=threshold,1,'first');

Ncycles=unique(idx_dist);
mov_exec_all=nan(Ndir*Npos*numel(Ncycles),1);
Nmax = size(scores,1);
Max_eigen_all=nan(Nmax,1);
idx_dir_new_all=nan(Nmax,1);
idx_pos_new_all=nan(Nmax,1);
idx_cycles_all=nan(Nmax,1);

counter = 1;

% for all conditions

for i_dist=1:numel(Ncycles)

    t_upto=sum(idx_dir==1 & idx_pos==1 & idx_dist==Ncycles(i_dist))/ms;
    mov_dur=(t_upto+timesmov(1)/ms-timesmov(2)/ms)*ms/Ncycles(i_dist); %average cycle duration

    this_cond = idx_dist==Ncycles(i_dist);
    % LDS for this condition
    [Max_eigen,idx_dir_new,idx_pos_new]=fit_LDS_PCs(scores(this_cond,:), ...
        idx_dir(this_cond),idx_pos(this_cond),explained,threshold,ndim^2*20);

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

%% Find dynamical epochs and plot them 
Detect_dynamics_epochs([real(Max_eigen_all),imag(Max_eigen_all)],idx_dir_new_all,idx_pos_new_all,idx_cycles_all,timesmov(1),mov_exec_all,do_plot);

end