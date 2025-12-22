function plot_fancy_errorbars(x,y,colours)
% PLOT_FANCY_ERRORBARS Plots jittered samples on the x axis with mean 
% and standard deviation. A default jitter of 0.1 is used.
%
%   PLOT_FANCY_ERRORBARS(x, y, colours)
%
%   INPUTS
%   ------
%   x        : [1 x Ncond] vector
%       X-axis positions corresponding to each condition.
%
%   y        : [Nsamples x Ncond] matrix
%       Sample values for each condition. Each column corresponds to one
%       condition, and each row to one sample. NaN values are ignored when
%       computing summary statistics.
%
%   colours  : [Ncond x 3] matrix
%       RGB colour values (in the range [0, 1]) used to plot the mean and
%       standard deviation for each condition.
%
%   OUTPUTS
%   -------
%   None
%       The function produces a plot in the current axes.
%
% Andrea Colins Rodriguez
% 21/12/2025

Ncond=numel(x);

Nsamples=size(y,1);

% jitter randomly so points are easily visible. 
if numel(x)>1
jitter=min(diff(x))*0.1;
else
    jitter=0.1;
end

% Plot on the same panel that is currently selected outside this function
hold on
for i_cond=1:Ncond

    plot(x(i_cond)+(rand(Nsamples,1)-0.5)*jitter,y(:,i_cond),'.','Color',[0.5 0.5 0.5])
    errorbar(x(i_cond),mean(y(:,i_cond),'omitnan'),std(y(:,i_cond),[],'omitnan'),'.','Color',colours(i_cond,:),'MarkerSize',12,'CapSize',0)

end

end