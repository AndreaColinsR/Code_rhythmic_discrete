function r=CCA_RNN_M1(scores_RNN,cond_idx_RNN,scores_M1,cond_idx_M1)
ndim=min([size(scores_RNN,2) 10]);
% interpolate RNN because of resizing

scores_RNN=scores_RNN(:,1:ndim);

cond_summ=unique(cond_idx_M1,'stable','rows');
%% match indexes of RNN to M1
new_scores_RNN=nan(size(scores_M1,1),ndim);

for i=1:size(cond_summ,1)
    this_new_cond_M1 = cond_idx_M1(:,1)==cond_summ(i,1) & cond_idx_M1(:,2)==cond_summ(i,2) & cond_idx_M1(:,3)==cond_summ(i,3);

    this_cond_RNN = cond_idx_RNN(:,1)==cond_summ(i,1) & cond_idx_RNN(:,2)==cond_summ(i,2) & cond_idx_RNN(:,3)==cond_summ(i,3);

    % can't interpolate across conditions, it needs to be done for each condition
    new_scores_RNN(this_new_cond_M1,:)=interp1(linspace(0,1,sum(this_cond_RNN))',scores_RNN(this_cond_RNN,:),linspace(0,1,sum(this_new_cond_M1))');
end

baseline=mean(new_scores_RNN(1:50,:));
new_scores_RNN=new_scores_RNN-baseline;



[~,~,r,RNN_CC,M1_CC] = canoncorr(new_scores_RNN,scores_M1(:,1:ndim));


do_plot=0;

if do_plot
%% check that they are still temporally aligned
figure


subplot(2,2,1)
plot(scores_M1(this_new_cond_M1,1))
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
