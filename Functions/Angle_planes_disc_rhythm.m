function Angle_dir = Angle_planes_disc_rhythm(states,idx_cycle,idx_Ncycle,idx_dir,idx_pos,do_plot)
% ANGLE_PLANES_DISC_RHYTHM Computes the angle between rotational planes of
% neural trajectories of rhythmic and discrete movements. Rotational planes are computed using jPCA 
% for rhythmic and discrete movement conditions. The output angle is the average across all movement directions and positions.
%
%
%   Angle_dir = ANGLE_PLANES_DISC_RHYTHM(states,idx_cycle,idx_Ncycle,...
%                                       idx_dir,idx_pos,do_plot)
%
%
%
%   INPUTS
%   ------
%   states : [T x N] matrix
%       Neural activity matrix (RNN or neural recording), where T is the number of time samples
%       and N is the number of units or neurons.
%
%   idx_cycle : [T x 1] vector
%       Movement cycle index. Rhythmic trials are identified as
%       idx_cycle == 7; discrete trials as idx_cycle < 1.
%
%   idx_Ncycle : [T x 1] vector
%       Cycle count index within a movement.
%
%   idx_dir : [T x 1] vector
%       Movement direction index for each movement.
%
%   idx_pos : [T x 1] vector
%       Movement position index for each movement.
%
%   do_plot : logical scalar
%       If do_plot equals 1, a 3-D PCA visualization of discrete and rhythmic trajectories
%       is generated. 
%
%   OUTPUT
%   ------
%   Angle_dir : scalar
%       Mean angle (in degrees) between the rhythmic and discrete rotational
%       planes, averaged across all directions and positions.
%
% Andrea Colins
% 17/12/2025



Ndir=max(idx_dir);
Npos=max(idx_pos);

Angle_dir=nan(Ndir,Npos);

for i_dir=1:Ndir
    for i_pos=1:Npos
        
        % Identify plane of rotation of rhythmic condition (during execution)
        Rhythmic=idx_cycle==7 & idx_Ncycle>1 & idx_Ncycle<7 & idx_dir==i_dir & idx_pos==i_pos;
        NelR=sum(Rhythmic);
        times_exec=1:NelR;
        [~,summary]=jPCA_prep_exec(states(Rhythmic,:),ones(NelR,1),ones(NelR,1),times_exec);
        
        % Identify plane of rotation of discrete condition (during execution)
        Discrete=idx_cycle<1 & idx_Ncycle==1 & idx_dir==i_dir & idx_pos==i_pos;
        NelD=sum(Discrete);
        times_exec=1:NelD;
        [~,summary2]=jPCA_prep_exec(states(Discrete,:),ones(NelD,1),ones(NelD,1),times_exec);

        % Compute angle between both planes of rotations
        Angle_dir(i_dir,i_pos) = subspace(summary.jPCs_highD(:,1:2),summary2.jPCs_highD(:,1:2))*180/pi;

        if do_plot && i_dir==1 && i_pos==1
            colour_cycle=plasma(5);

            Discrete=idx_cycle<1 & idx_dir==i_dir & idx_pos==i_pos;
            NelD=sum(Discrete);

            [~,scores]=pca([states(Discrete,:);states(Rhythmic,:)]);
            plot3(scores(1:NelD,1),scores(1:NelD,2),scores(1:NelD,3),'Color',colour_cycle(1,:))

            hold on

            exec_disc=find(idx_Ncycle(Discrete)>0);
            plot3(scores(exec_disc,1),scores(exec_disc,2),scores(exec_disc,3),'Color',colour_cycle(1,:),'LineWidth',2)
            plot3(scores(1,1),scores(1,2),scores(1,3),'ok')

            plot3(scores(NelD+1:end,1),scores(NelD+1:end,2),scores(NelD+1:end,3),'Color',colour_cycle(5,:),'LineWidth',2)
            title(['Angle = ' num2str(Angle_dir(i_dir,i_pos),'%.2f')])

        end

    end
end

Angle_dir=mean(Angle_dir(:));
end