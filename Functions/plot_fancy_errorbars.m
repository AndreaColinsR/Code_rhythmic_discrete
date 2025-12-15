function plot_fancy_errorbars(x,y,colours)

Ncond=numel(x);

Nsamples=size(y,1);

if numel(x)>1
jitter=min(diff(x))*0.1;
else
    jitter=0.1;
end
hold on
for i_cond=1:Ncond

    plot(x(i_cond)+(rand(Nsamples,1)-0.5)*jitter,y(:,i_cond),'.','Color',[0.5 0.5 0.5])
    errorbar(x(i_cond),mean(y(:,i_cond),'omitnan'),std(y(:,i_cond),[],'omitnan'),'.','Color',colours(i_cond,:),'MarkerSize',12,'CapSize',0)

end

end