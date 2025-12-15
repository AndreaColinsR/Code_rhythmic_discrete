function plot_trajectories_colour(ax,trajectory,nbins)
colours=pink(2*nbins);
colours=colours(1:nbins,:);
segments=round(size(trajectory,1)./nbins);
plot3(ax,trajectory(1,1),trajectory(1,2),trajectory(1,3),'o','MarkerFaceColor',colours(1,:),'MarkerEdgeColor',colours(1,:))
hold on
for i=1:nbins-1
    idx=(i-1)*segments+1:i*segments+1;
    plot3(ax,trajectory(idx,1),trajectory(idx,2),trajectory(idx,3),'Color',colours(i,:))
end
plot3(ax,trajectory(idx(end):end,1),trajectory(idx(end):end,2),trajectory(idx(end):end,3),'Color',colours(end,:))
end