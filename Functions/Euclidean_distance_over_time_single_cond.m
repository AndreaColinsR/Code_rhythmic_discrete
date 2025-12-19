function [total_dist,comb]=Euclidean_distance_over_time_single_cond(score,idx_condition)
Ncond=max(idx_condition);
Nsamples=nan(Ncond,1);

for i_cond=1:Ncond

    Nsamples(i_cond)=sum(idx_condition==i_cond);
end

%cut all trajectories to the minimum number of samples

Nsamples=min(Nsamples);
Ncomb=sum(1:Ncond-1);
total_dist=nan(Nsamples,Ncomb);
comb=nan(Ncomb,2);
counter=1;
for i_cond=1:Ncond

    idx1=find(idx_condition==i_cond);
    % subplot(4,4,7)
    % hold on
    % plot3(score(idx1(1:Nsamples),1),score(idx1(1:Nsamples),2),score(idx1(1:Nsamples),3))

    for j_cond=i_cond+1:Ncond
        idx2=find(idx_condition==j_cond);
        
        % subplot(4,4,7)
        % plot3(score(idx2(1:Nsamples),1),score(idx2(1:Nsamples),2),score(idx2(1:Nsamples),3))
        % hold off

        total_dist(:,counter)=sqrt(sum((score(idx1(1:Nsamples),:)-score(idx2(1:Nsamples),:)).^2,2));
        comb(counter,:)=[i_cond,j_cond];
        counter=counter+1;
      
    end
end
end