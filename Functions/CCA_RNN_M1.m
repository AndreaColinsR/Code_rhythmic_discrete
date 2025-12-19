function r = CCA_RNN_M1(scores_RNN,cond_idx_RNN,scores_Neural,cond_idx_Neural)
% CCA_RNN_M1 Performs canonical correlation analysis (CCA) between low-dimensional
% representations of RNN activity and neural activity.
%
%   CCA_RNN_M1(scores_RNN,cond_idx_RNN,scores_Neural,cond_idx_Neural)
%
%   INPUTS
%   ------
%
%   scores_RNN  : [T_RNN x D] matrix
%       Low-dimensional RNN activity (e.g. PCA scores), where T_RNN is the
%       number of RNN time samples and D is the number of dimensions.
%
%   cond_idx_RNN : [T_RNN x 3] matrix
%       Condition index for RNN activity. Columns correspond to task
%       variables (e.g. direction, position and distance).
%
%   scores_Neural   : [T_Neural x D] matrix
%       Low-dimensional M1 neural activity (e.g. PCA scores), where T_Neural is
%       the number of cortical time samples and D is the number of
%       dimensions.
%
%   cond_idx_Neural : [T_Neural x 3] matrix
%       Condition index for M1 activity. Columns correspond to task
%       variables (e.g. direction, position and distance).
%
%   OUTPUTS
%   -------
%
%   r           : [N x 1] vector
%       Canonical correlation coefficients between RNN and M1 activity,
%       ordered from highest to lowest correlation.
%
%
% Andrea Colins Rodriguez
% 19/12/2025

ndim=min([size(scores_RNN,2) 10]);
scores_RNN=scores_RNN(:,1:ndim);

cond_summ=unique(cond_idx_Neural,'stable','rows');

%% Match indexes of RNN to M1
new_scores_RNN=nan(size(scores_Neural,1),ndim);

for i=1:size(cond_summ,1)

    this_new_cond_M1 = cond_idx_Neural(:,1)==cond_summ(i,1) & cond_idx_Neural(:,2)==cond_summ(i,2) & cond_idx_Neural(:,3)==cond_summ(i,3);

    this_cond_RNN = cond_idx_RNN(:,1)==cond_summ(i,1) & cond_idx_RNN(:,2)==cond_summ(i,2) & cond_idx_RNN(:,3)==cond_summ(i,3);
    
    % Interpolate RNN so that its neural activity (dt = 10) has the same temporal resolution
    % that cortical activity (dt = 1).
    % can't interpolate across conditions, it needs to be done for each condition
    new_scores_RNN(this_new_cond_M1,:)=interp1(linspace(0,1,sum(this_cond_RNN))',scores_RNN(this_cond_RNN,:),linspace(0,1,sum(this_new_cond_M1))');
end

% we define the baseline of the RNN activity as the average in the first
% 50 ms
baseline=mean(new_scores_RNN(1:50,:));
new_scores_RNN=new_scores_RNN-baseline;

[~,~,r] = canoncorr(new_scores_RNN,scores_Neural(:,1:ndim));


% Sanity check
do_plot=0;

if do_plot
%% check that they are still temporally aligned
figure
[~,~,~,RNN_CC,M1_CC] = canoncorr(new_scores_RNN,scores_Neural(:,1:ndim));

subplot(2,2,1)
plot(scores_Neural(this_new_cond_M1,1))
hold on
plot(new_scores_RNN(this_new_cond_M1,1))
legend('Neural','RNN')
title('PC1 befor aligment')

subplot(2,2,2)
plot(1:ndim,r)

title('CCs')

subplot(2,2,3)
hold on
plot(M1_CC(this_new_cond_M1,1))
plot(RNN_CC(this_new_cond_M1,1))
legend('Neural','RNN')
title('CC1 after aligment')


subplot(2,2,4)
hold on 
plot3(M1_CC(this_new_cond_M1,1),M1_CC(this_new_cond_M1,2),M1_CC(this_new_cond_M1,3))
plot3(RNN_CC(this_new_cond_M1,1),RNN_CC(this_new_cond_M1,2),RNN_CC(this_new_cond_M1,3))
legend('Neural','RNN')
title('CCs after aligment')
end


end
