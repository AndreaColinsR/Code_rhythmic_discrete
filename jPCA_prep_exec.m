function [scorestruct,summary,Data]=jPCA_prep_exec(FR,idx_pos,idx_dir,selected_times)
%% aim: to test if the planes of rotation during prep and execution are the same

Ndim=4;
Npos=max(idx_pos);
Ndir=max(idx_dir);

%% for jpca
jPCA_params.params=false;
jPCA_params.numPCs = Ndim;  % default anyway, but best to be specific
jPCA_params.suppressBWrosettes = true;  % these are useful sanity plots, but lets ignore them for now
jPCA_params.suppressHistograms = true;  % these are useful sanity plots, but lets ignore them for now
jPCA_params.suppressText=true;
jPCA_params.softenNorm=5;
jPCA_params.normalize=false;


%total_time=length(idx);


i_cond=1;
for i_pos=1:Npos
    for i_dir=1:Ndir
        idx=find(idx_pos==i_pos & idx_dir==i_dir);
        Data(i_cond).A=FR(idx(selected_times),:);
        Data(i_cond).times=(selected_times)';
        i_cond=i_cond+1;
    end
end
%Ncond=i_cond-1;

[scorestruct,summary] = jPCA(Data,selected_times, jPCA_params);
[summary.varCaptEachJPC,idxtmp]=sort(summary.varCaptEachJPC,'descend');
summary.jPCs=summary.jPCs(:,idxtmp);
summary.jPCs_highD=summary.jPCs_highD(:,idxtmp);
% figure
% hold on
% for i_cond=1:Ncond
% jpcs=scorestruct(i_cond).tradPCAproj(:,1:Ndim);
% 
% plot3(jpcs(:,1),jpcs(:,2),jpcs(:,3))
% 
% plot3(jpcs(1,1),jpcs(1,2),jpcs(1,3),'o')
% end


%tmpdim=find(cumsum(summary.varCaptEachJPC)>=(varExpTh/100),1,'First');
end