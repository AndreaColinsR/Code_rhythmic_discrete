function plot_std(x,y,std_dev,Colour)
% PLOT_STD Plots a curve with a shaded standard deviation envelope (± one standard deviation).
% Any NaN values in the input data are excluded prior to plotting.
%
%   PLOT_STD(x, y, std_dev, Colour)
%
%   INPUTS
%   ------
%   x        : [1 x N] vector
%       X-axis values corresponding to the mean data points.
%
%   y        : [1 x N] vector
%       Mean values to be plotted as a continuous curve.
%
%   std_dev  : [1 x N] vector
%       Standard deviation associated with each mean value.
%
%   Colour   : [1 x 3] vector
%       RGB colour specification (values in the range [0, 1]) used for both
%       the shaded region and the mean curve.
%
%   OUTPUTS
%   -------
%   None
%       The function produces a plot in the current axes.
%
%
% Andrea Colins
% 21/12/2025



select=~isnan(y);
x=x(select);
std_dev=std_dev(select);
y=y(select);

% check x and y are a row
x=x(:)';
y=y(:)';
std_dev=std_dev(:)';

curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

%% plot shaded area 
fill(x2, inBetween,Colour,'FaceAlpha',.5,'EdgeColor',Colour,'EdgeAlpha',.5)
hold on
%% plot curve on top (mean)
plot(x, y, 'Color',Colour, 'LineWidth', 2);
end