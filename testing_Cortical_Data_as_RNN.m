function testing_Cortical_Data_as_RNN(region_name,figW,plot_supp_figs)

animal={'Cousteau','Drake'}; % 'E' or 'F'
dists=[0.5 1 2 4 7];
Ndists=numel(dists);
Nanimals=numel(animal);
colour_dist=plasma(Ndists);
Nplots = 4;

column = 3;

Mean_over_time=nan(1,Nanimals);
Mean_dist_An=nan(6000,Ndists,Nanimals);
Init_cond_t_all=nan(8,Ndists);

if plot_supp_figs.do_plot == 1
    if strcmp(region_name,'M1')
        prep_fig = plot_supp_figs.Prep_M1;
        dpca_fig = plot_supp_figs.dPCA_M1;
    else
        prep_fig = plot_supp_figs.Prep_SMA;
        dpca_fig = plot_supp_figs.dPCA_SMA;
    end
else
    prep_fig = [];
    dpca_fig = [];
end

for i_animal=1:Nanimals

    load(['.\Output_files\scores_' animal{i_animal} '_'  region_name '.mat'],'FR','scores','idx_dir','idx_pos','idx_dist','exec','idx_Ncycle')
     
    %% PLot example trajectories
    if (strcmp(animal{i_animal},'Drake') && strcmp(region_name,'M1')) || (strcmp(animal{i_animal},'Cousteau') && strcmp(region_name,'SMA'))
        figure(figW)
        i_dir = 2;
        i_pos = 1;
        idx_rhythmic=idx_dir==i_dir & idx_pos==i_pos & idx_dist==7;
        idx_discrete=idx_dir==i_dir & idx_pos==i_pos & idx_dist==0.5;
        scores_rhythmic_M1=scores(idx_rhythmic,:);
        scores_discrete_M1=scores(idx_discrete,:);
        view_angle=[48,21];

        % plotting Neural trajectories
        subplot(Nplots,Nplots,7)
        hold on
        plot3(scores_rhythmic_M1(:,1),scores_rhythmic_M1(:,2),scores_rhythmic_M1(:,3),'Color',colour_dist(end,:))
        plot3(scores_discrete_M1(:,1),scores_discrete_M1(:,2),scores_discrete_M1(:,3),'Color',colour_dist(1,:))
        plot3(scores_rhythmic_M1(1,1),scores_rhythmic_M1(1,2),scores_rhythmic_M1(1,3),'o','MarkerFaceColor',colour_dist(end,:),'MarkerEdgeColor',colour_dist(end,:))
        plot3(scores_discrete_M1(1,1),scores_discrete_M1(1,2),scores_discrete_M1(1,3),'o','MarkerFaceColor',colour_dist(1,:),'MarkerEdgeColor',colour_dist(1,:))

        plot3(scores_rhythmic_M1(1000:end-400,1),scores_rhythmic_M1(1000:end-400,2),scores_rhythmic_M1(1000:end-400,3),'Color',colour_dist(end,:),'LineWidth',2)
        plot3(scores_discrete_M1(1000:end-400,1),scores_discrete_M1(1000:end-400,2),scores_discrete_M1(1000:end-400,3),'Color',colour_dist(1,:),'LineWidth',2)

        view(view_angle)
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
    end

    %% Plot preparation (supp)
    if plot_supp_figs.do_plot == 1
        if (strcmp(animal{i_animal},'Cousteau') && strcmp(region_name,'M1')) || (strcmp(animal{i_animal},'Drake') && strcmp(region_name,'SMA')) 
                figure(prep_fig)
            plot_preparation_all_conditions(FR,idx_dir, idx_pos, idx_dist,exec,column)
        end
    end

    %% Plot initial conditions in SMA
    if strcmp(region_name,'SMA') && strcmp(animal{i_animal},'Cousteau')
        figure(figW)
        plot_init = 1;
        subplot(4,4,8+column)
    else
        plot_init = 0;
    end

    %% Plot Linear Dynamical Sytem (LDS) and angle results
    if strcmp(region_name,'M1') && strcmp(animal{i_animal},'Cousteau')
        plot_LDS =1;
    else
        plot_LDS = 0;
    end

    [Angle_disc_rhythm,Init_cond_t,Dist2Att] = RNNs_predictions(FR,idx_dir,idx_pos,idx_dist,exec,idx_Ncycle,column,plot_init,plot_LDS,plot_supp_figs.LDS);

    Init_cond_t_all(4*(i_animal-1)+1:4*i_animal,:) = Init_cond_t;


    [~,prepdata]=get_prep_exec_after_FR(FR,idx_pos,idx_dir,idx_dist,idx_Ncycle);

    %% dPCA
    do_plot = 0;
    if plot_supp_figs.do_plot == 1
        if (strcmp(animal{i_animal},'Drake') && strcmp(region_name,'M1')) || (strcmp(animal{i_animal},'Drake') && strcmp(region_name,'SMA'))
            do_plot = 1;
            figure(dpca_fig)
        else
            do_plot = 0;
        end
    end

    [~,Per_prep] = dPCA_across_conditions(prepdata.FR,prepdata.ndir,prepdata.npos,prepdata.ndist,do_plot,column);

    if plot_supp_figs.do_plot == 1
        figure(dpca_fig)
        subplot(5,3,column*5)
        hold on
        plot(1:4+(i_animal-1)*0.2,Per_prep,'.-b')
        xticks(1:4)
        xticklabels({'N cycle','Direction','Pos','CI'})
        figure(figW)

    end


    Mean_dist_An(1:size(Dist2Att,1),:,i_animal)=squeeze(mean(Dist2Att(:,:,:),2));
    Mean_over_time(i_animal)=min(mean(Dist2Att(:,:,1),2));
    MinDist=squeeze(min(Dist2Att(:,:,1:2)));

    %% Ploting distances to between trajectories of different conditions

    ndims=6;
    Dist_all=compare_traj_directions(prepdata.scores(:,1:ndims),prepdata.ndir,prepdata.npos,prepdata.ndist);
    prept=size(Dist_all,1);
    t=-prept:-1;

    subplot(Nplots,Nplots,15)
    hold on
    plot(t,Dist_all(:,1),'Color',[0 0 0])
    plot(t,Dist_all(:,2),'Color',[0.5 0.6 0.9])
    plot(t,Dist_all(:,3),'Color',[0.5 0.5 0.5])


    save(['.\Output_files\scores_' animal{i_animal} '_'  region_name '.mat'],'Per_prep','Angle_disc_rhythm','Init_cond_t','MinDist','Mean_over_time','Dist_all','-append')
end

%% Plotting distance to the attractor (M1)
if strcmp(region_name,'M1')
    subplot(Nplots,Nplots,11)
    hold on

    t2=-prept:size(Mean_dist_An)-prept-1;
    plot([t2(1) t2(end)],[mean(Mean_over_time) mean(Mean_over_time)],'Color',[0.5 0.5 0.5])

    for i_dist=1:numel(dists)
        plot(t2,mean(Mean_dist_An(:,i_dist),3),'Color',colour_dist(i_dist,:))
    end

else
    subplot(Nplots,Nplots,12)
    hold on
    errorbar(mean(Init_cond_t_all,'omitnan'),[0.5 1 2 4 7],std(Init_cond_t_all,'omitnan'),'.-k','horizontal')

end

%% LDS (Supp M1)
if strcmp(region_name,'M1') && plot_supp_figs.do_plot == 1
    figure(plot_supp_figs.LDS)
    subplot(2,3,column)
    Main_LDS_Cortical('Cousteau',1)

end