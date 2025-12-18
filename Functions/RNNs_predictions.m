function [Angle_disc_rhythm,Init_cond_t,Dist2Att] = RNNs_predictions(states,idx_dir,idx_pos,idx_cycle,exec,idx_Ncycle,do_plot_pred,column,plot_init,fig_supp,plot_angle,fig_angle)
%%
% Angle_dir angle between rotational planes of discrete and rhythmic
% movements
% Dist2Att is the distance to the limit cycle attractor

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

%% if plotting an example of neural activity during preparation

if do_plot_pred

Ndir=2;
Npos=2;

NumberCyles=unique(idx_cycle);
Ncycle=numel(NumberCyles);

colour_dir=[[115, 80, 185];[87, 204, 153]]./255;%[[194 165 207];[166 219 160]]./255;
colour_ndist=plasma(Ncycle);

    figure(fig_supp)
    [Ntimes,Nunits] = size(states);
    counter = 1;

    idxCycle_prep=nan(Ntimes,1);
    idxDir_prep=nan(Ntimes,1);
    idxPos_prep=nan(Ntimes,1);
    FRicycle_prep_all=nan(Ntimes,Nunits);

    for i_pos=1:Npos
        for i_dir=1:Ndir


            for i_cycle=1:Ncycle

                this_cond=find(idx_dir==i_dir & idx_pos==i_pos & idx_cycle==NumberCyles(i_cycle));
                this_prep=1:find(exec(this_cond)==1,1,'first');

                FRicycle=states(this_cond(this_prep),:);

                Nrows = numel(this_cond(this_prep));

                idxCycle_prep(counter:counter+Nrows-1,:)=ones(Nrows,1)*i_cycle;
                idxDir_prep(counter:counter+Nrows-1,:)=ones(Nrows,1)*i_dir;
                idxPos_prep(counter:counter+Nrows-1,:)=ones(Nrows,1)*i_pos;
                FRicycle_prep_all(counter:counter+Nrows-1,:)=FRicycle;

                counter=counter+Nrows;
            end


        end
    end

    idx_nan = isnan(idxCycle_prep);
    idxCycle_prep(idx_nan,:) = [];
    idxDir_prep(idx_nan,:) = [];
    idxPos_prep(idx_nan,:) = [];
    FRicycle_prep_all(idx_nan,:) = [];

    [~,PCAPrep]=pca(FRicycle_prep_all);

    for i_pos=1:Npos
        for i_dir=1:Ndir
            for i_cycle=1:Ncycle

                this_cond=find(idxCycle_prep==i_cycle & idxDir_prep==i_dir & idxPos_prep==i_pos);
                subplot(2,3,3+column)
                hold on
                plot3(PCAPrep(this_cond(this_prep),1),PCAPrep(this_cond(this_prep),2),PCAPrep(this_cond(this_prep),3),'Color',colour_ndist(i_cycle,:),'LineWidth',5/i_cycle)
                plot3(PCAPrep(this_cond(this_prep(1)),1),PCAPrep(this_cond(this_prep(1)),2),PCAPrep(this_cond(this_prep(1)),3),'o','Color',colour_ndist(i_cycle,:))

                subplot(2,3,column)
                hold on
                plot3(PCAPrep(this_cond(this_prep),1),PCAPrep(this_cond(this_prep),2),PCAPrep(this_cond(this_prep),3),'Color',colour_dir(i_dir,:)./sqrt(i_pos),'LineWidth',5/i_cycle)
                plot3(PCAPrep(this_cond(this_prep(1)),1),PCAPrep(this_cond(this_prep(1)),2),PCAPrep(this_cond(this_prep(1)),3),'o','Color',colour_dir(i_dir,:)./sqrt(i_pos))


            end
        end
    end

end


end