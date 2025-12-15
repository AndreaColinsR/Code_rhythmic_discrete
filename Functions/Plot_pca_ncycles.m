function Plot_pca_ncycles(scores,idx_dir,idx_pos,idx_dist,varargin)
if nargin<5
    figure
    ax=[subplot(1,3,1),subplot(1,3,2),subplot(1,3,3)];
else
    ax=varargin{1};
end

dt=10;
if size(scores,1)<10000
    %RNN
    dt=1;
end

%Ndir=max(idx_dir);
%Npos=max(idx_pos);
Ndist=[min(idx_dist) max(idx_dist)];
NNdist=numel(Ndist);
colour_dist=plasma(NNdist);

for i_dist=1:NNdist
    idx=find(idx_dir==2 & idx_pos==1 & idx_dist==Ndist(i_dist));

    %% over time
    axes(ax(1+i_dist))
    hold on
    plot_trajectories_colour(ax(1+i_dist),scores(idx,:),round(numel(idx)/40))
    


%     ax2=subplot(2,NNdist,NNdist+1);
%     hold on
%     plot_trajectories_colour(ax2,scores(idx,:),round(numel(idx)/40))

%     subplot(2,NNdist,i_dist)
%     hold on
%     plot3(scores(idx,1),scores(idx,2),scores(idx,3),'Color',colour_dist(i_dist,:))
%     plot3(scores(idx(1),1),scores(idx(1),2),scores(idx(1),3),'.','MarkerSize',12,'Color',colour_dist(i_dist,:))
%     plot3(scores(idx(100*dt:end-40*dt),1),scores(idx(100*dt:end-40*dt),2),scores(idx(100*dt:end-40*dt),3),'Color',colour_dist(i_dist,:),'LineWidth',3)
%     view(64.5239,15.9280)


    % coloured by cycle
    axes(ax(1));
    hold on
    plot3(scores(idx,1),scores(idx,2),scores(idx,3),'Color',colour_dist(i_dist,:))
    plot3(scores(idx(1),1),scores(idx(1),2),scores(idx(1),3),'.','MarkerSize',30,'Color',colour_dist(i_dist,:))
    plot3(scores(idx(100*dt:end-40*dt),1),scores(idx(100*dt:end-40*dt),2),scores(idx(100*dt:end-40*dt),3),'Color',colour_dist(i_dist,:),'LineWidth',3)
    
    view(64.5239,15.9280)

   
    %% tmp to find example 




%     % monkey 2
%     %view(162.6198,3.4843)
end

end