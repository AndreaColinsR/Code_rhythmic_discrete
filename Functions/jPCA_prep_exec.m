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

template_struct.A = [];
template_struct.times = [];
Data = repmat(template_struct, 1,Ndir*Npos);

i_cond=1;
for i_pos=1:Npos
    for i_dir=1:Ndir
        idx=find(idx_pos==i_pos & idx_dir==i_dir);
        Data(i_cond).A=FR(idx(selected_times),:);
        Data(i_cond).times=(selected_times)';
        i_cond=i_cond+1;
    end
end


[scorestruct,summary] = jPCA_new_version(Data,selected_times, jPCA_params);
[summary.varCaptEachJPC,idxtmp]=sort(summary.varCaptEachJPC,'descend');
summary.jPCs=summary.jPCs(:,idxtmp);
summary.jPCs_highD=summary.jPCs_highD(:,idxtmp);

end