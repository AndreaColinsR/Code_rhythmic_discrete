function [Angle_disc_rhythm,Init_cond_t,Dist2Att] = RNNs_predictions(states,idx_dir,idx_pos,idx_cycle,exec,idx_Ncycle,column,plot_init,plot_angle,fig_angle)
%  RNN_PREDICTIONS Computes three predictions of the dynamics of RNNs or
%  cortical activity.
%
%
%
%   [Angle_disc_rhythm,Init_cond_t,Dist2Att] = RNN_PREDICTIONS(states,idx_dir,idx_pos,...
%                   idx_cycle,exec,idx_Ncycle,column,plot_init,plot_angle,fig_angle)
%
%
%
%   INPUTS
%   ------
%   states : [T x N] matrix
%       Neural activity matrix (RNN or neural recording), where T is the number of time samples
%       and N is the number of units or neurons.
%
%   idx_dir     : [T x 1] vector of movement or task direction labels.
%
%   idx_pos     : [T x 1] vector of position or task condition labels.
%
%   idx_cycle   : [T x 1] vector indicating the total distance (Number of cycles) the animal
%   covers in the trial [0.5 1 2 4 7]
%
%   exec        : [T x 1] binary vector indicating if animal is performing
%   a movement at time t (1) or not (0). 
%
%   idx_Ncycle  : [T x 1] vector indicating the current cycle number performed within a
%                 movement.
%
%   column      : scalar
%                 indicates de column that the plot will be placed (1 (Same-control RNN),2 (Diff-control RNN) or 3 (cortical))
%
%   plot_init   : logical scalar
%                 If do_plot equals 1, a 3-D PCA visualization of discrete and rhythmic neural trajectories
%                 is generated. This plots shows the initial conditions at movement
%                 execution (SMA)
%
%   plot_angle  : logical scalar
%                 If do_plot equals 1, a 3-D PCA visualization of discrete and rhythmic neural trajectories
%                 is generated. This plots highlights the difference in planes of
%                 rotations between neural trajectories.
%
%   fig_angle : figure name
%               if plot_angle equals 1, the plots will be places in this
%               figure in the column indicated by the variable column. The
%               figure will be empty otherwise. 
%       . 
%
%   OUTPUT
%   ------
%   Angle_dir : scalar
%       Mean angle (in degrees) between the rhythmic and discrete rotational
%       planes, averaged across all directions and positions.
%
%   Dist2Att  : [Ntime x Ncond x Ncycle] array containing the normalized
%                 distance to the attractor.
%                 - Ntime  : maximum number of time points (with buffer)
%                 - Ncond  : direction × position conditions
%                 - Ncycle : unique cycle types
%
% Andrea Colins Rodriguez
% 18/12/2025

% soft normalise states
states=states-mean(states);
states=states./repmat(range(states)+5,size(states,1),1);

%% computing the 3 outputs:
% 1. Distance to limit cycle of rhythmic movements (for M1)
Dist2Att = distance_to_attractor(states,idx_cycle,idx_Ncycle,idx_dir,idx_pos);

% 2. Position of the neural state at mov onset (initial condition) relative
% to the helix (from SMA)
Init_cond_t = initial_cond_helix(states,idx_dir,idx_pos,idx_cycle,idx_Ncycle,exec,plot_init);

% 3. Angles between planes of rotations of neural activity during the execution of
% and rhythmic movements (M1)
if plot_angle
    figure(fig_angle)
    subplot(2,4,4+column)
end

Angle_disc_rhythm = Angle_planes_disc_rhythm(states,idx_cycle,idx_Ncycle,idx_dir,idx_pos,plot_angle);

end