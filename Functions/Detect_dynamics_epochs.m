function [epochs,epochs_all]=Detect_dynamics_epochs(eigs_traj,idx_dir,idx_speed,idx_cycles,tstart,t_end,do_plot)
% Detect_dynamics_epochs Identifies dynamical epochs from LDS eigenvalue
% trajectories.
%
%   [epochs,epochs_all] = Detect_dynamics_epochs(eigs_traj,idx_dir, ...
%       idx_speed,idx_cycles,tstart,t_end,do_plot)
%
%   This function classifies time-varying dynamical regimes based on the
%   eigenvalues of locally fitted linear dynamical systems (LDS). Each time
%   point is assigned to a discrete dynamical epoch reflecting expansion,
%   contraction and/or rotational dynamics. Epochs are aggregated across
%   behavioural conditions and can be visualised relative to movement
%   timing. Rotational dynamics are defined by a minimum imaginary component
%   corresponding to a maximum oscillatory period of 10 s.Expansion and 
%   contraction thresholds are defined on the real part of the eigenvalues.
%
%   INPUTS
%   ------
%
%   eigs_traj : [T x 2] matrix
%       Trajectory of dominant LDS eigenvalues, where the first column
%       contains the real part and the second column contains the imaginary
%       part of the eigenvalues.
%
%   idx_dir : [T x 1] vector
%       Direction index for each time point.
%
%   idx_speed : [T x 1] vector
%       Speed or position index for each time point.
%
%   idx_cycles : [T x 1] vector
%       Cycle index for each time point, defining repeated movement cycles.
%
%   tstart : scalar
%       Start time (in milliseconds) relative to movement onset.
%
%   t_end : [Ncond x 1] vector
%       End time (in milliseconds) of movement execution for each condition.
%
%   do_plot : logical
%       Flag indicating whether dynamical epochs should be visualised.
%
%   OUTPUTS
%   -------
%
%   epochs : [T x 1] vector
%       Dynamical epoch label assigned to each time point.
%
%   epochs_all : [L x Ncond] matrix
%       Dynamical epoch labels aligned in time and concatenated across all
%       conditions, where L is the maximum trajectory length.
%
%   DESCRIPTION
%   -----------
%
%   Dynamical regimes are defined using thresholds on the real and imaginary
%   parts of the dominant eigenvalues. Each time point is classified into
%   one of five categories:
%
%       1 – contraction only
%       2 – rotation with contraction
%       3 – rotation only
%       4 – rotation with expansion
%       5 – expansion only
%
%   Epoch labels are grouped by direction, speed and cycle condition to
%   produce condition-aligned dynamical trajectories. When requested, the
%   function visualises the resulting epochs using a colour-coded map
%   aligned to movement onset and offset.
%
%
% Andrea Colins Rodriguez
% 22/12/2025

ms=1000;
max_period=10*ms;
threshold_rot=2*pi/max_period;
thr_exp=0.0015;

Ndir=max(idx_dir);
Nspeed=max(idx_speed);
Ncycles=numel(unique(idx_cycles));

Ncond=Ndir*Nspeed*Ncycles;
% hard coding the maximum length
max_l=round(max(t_end)-tstart+1000);
epochs_all=zeros(max_l,Ncond);
epochs=zeros(size(eigs_traj,1),1);


expan=eigs_traj(:,1)>0 & (isnan(eigs_traj(:,2)) | eigs_traj(:,2)<threshold_rot);
contra=eigs_traj(:,1)<0 & (isnan(eigs_traj(:,2)) | eigs_traj(:,2)<threshold_rot);
rot=abs(eigs_traj(:,1))<thr_exp & eigs_traj(:,2)>=threshold_rot;
rot_exp=eigs_traj(:,1)>thr_exp & eigs_traj(:,2)>=threshold_rot;
rot_cont=eigs_traj(:,1)<-thr_exp & eigs_traj(:,2)>=threshold_rot;


epochs(contra)=1;
epochs(rot_cont)=2;
epochs(rot)=3;
epochs(rot_exp)=4;
epochs(expan)=5;


i_cond=1;
for i_cycles=1:Ncycles
    for i_dir=1:Ndir

        for i_speed=1:Nspeed

            idx=idx_dir==i_dir & idx_speed==i_speed & idx_cycles==i_cycles;
            epochs_all(1:sum(idx),i_cond)=epochs(idx);

            i_cond=i_cond+1;
        end
    end
end

if do_plot

colour_epoch=[255 255 255;0,90,116;74,148,159;196,196,196; 244,119,127;147,0,58]./255;

imagesc([tstart max_l+tstart],[1 Ncond],epochs_all')
hold on
plot([0 0],[0 Ncond],'Color',[0 0 0],'LineWidth',2)
t_2=reshape([t_end,t_end]',Ncond*2,1);
tendexec=reshape([(1:Ncond)-0.5;(1:Ncond)+0.5],Ncond*2,1);
plot(t_2,tendexec,'Color',[0 0 0],'LineWidth',2)
colormap(colour_epoch)
clim([0 5])
box off
xlabel('Time to movement onset [ms]')
ylabel('Condition')
cbh = colorbar ; %Create Colorbar
cbh.Ticks = linspace(0.5,4.5,6); %Create 8 ticks from zero to 1
cbh.TickLabels = {'','only contraction','rotation+ contraction','only rotation','rotation+ expansion','only expansion'} ; 
end
end
