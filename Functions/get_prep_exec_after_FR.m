function [execdata,prepdata]=get_prep_exec_after_FR(FR,idx_pos,idx_dir,idx_dist,idx_Ncycle)


Ndir=max(idx_dir);
Npos=max(idx_pos);
Ndist=unique(idx_dist);

prep=[];
exec=[];

ndist_prep=[];
ndir_prep=[];
npos_prep=[];



ncycle_exec=[];
ndist_exec=[];
ndir_exec=[];
npos_exec=[];
npos_endexec=[];

ndist_exec_all=[];
ndir_exec_all=[];
npos_exec_all=[];



tmp=600;

for i_dist=1:numel(Ndist)
    for i_dir=1:Ndir
        for i_pos=1:Npos
            %% for each condition select prep
            cond=find(idx_dir==i_dir & idx_pos==i_pos & idx_dist==Ndist(i_dist));

            idx_onset=find(idx_Ncycle(cond)>0,1,"first")-1;
            idx_onset2= idx_onset; %start of execution for dpca

            prep=[prep;FR(cond(1:idx_onset),:)];
            ndist_prep=[ndist_prep;ones(idx_onset,1)*i_dist];
            ndir_prep=[ndir_prep;ones(idx_onset,1)*i_dir];
            npos_prep=[npos_prep;ones(idx_onset,1)*i_pos];

            %% for each condition select exec
            idx_mov_end=find(~isnan(idx_Ncycle(cond)),1,'last');
            Nexec=numel(cond(idx_onset+1:idx_mov_end));
            ndist_exec_all=[ndist_exec_all;ones(Nexec,1)*i_dist];
            ndir_exec_all=[ndir_exec_all;ones(Nexec,1)*i_dir];
            npos_exec_all=[npos_exec_all;ones(Nexec,1)*i_pos];
            
            %% for each condition select exec of first half cycle for dpca
            % start of second cycle or end of sequence
            if Ndist(i_dist)<1
                % for half a cycle, take the entire movement
                idx_offset=idx_mov_end;
                    
               
            elseif Ndist(i_dist)==1
                 % for one cycle, take the first half cycle
                firstcycledur=find(isnan(idx_Ncycle(cond(idx_onset+1:end))),1,"first");
                idx_offset=idx_onset+round(firstcycledur/2);


            else
                % for more than one cycle, take the first half cycle
                firstcycledur=find(idx_Ncycle(cond(idx_onset+1:end))>1,1,"first");
                idx_offset=idx_onset+round(firstcycledur/2);
            end
            

            FRfirstHalf=interp1(linspace(0,1,idx_offset-idx_onset2+1)',FR(cond(idx_onset2:idx_offset),:),linspace(0,1,tmp)');

         

            exec=[exec;FRfirstHalf];

            ndist_exec=[ndist_exec;ones(tmp,1)*i_dist];
            ndir_exec=[ndir_exec;ones(tmp,1)*i_dir];
            npos_exec=[npos_exec;ones(tmp,1)*i_pos];
            
            % for end of exec, switch the init pos index
            if Ndist(i_dist)>=1
                npos_endexec=[npos_endexec;ones(tmp,1)*i_pos];
            else    
                if i_pos==1
                    npos_endexec=[npos_endexec;ones(tmp,1)*2];
                else
                    npos_endexec=[npos_endexec;ones(tmp,1)*1];
                end
            end

            ncycle_exec=[ncycle_exec;idx_Ncycle(cond)];
        end
    end
end
selected_neurons=sum(FR)>1e-6;
FR=FR(:,selected_neurons);


% soft normalise
norm_factor=range(FR)+5;
prep=prep(:,selected_neurons)./norm_factor;
exec=exec(:,selected_neurons)./norm_factor;

[coeff_prep,score_prep]=pca(prep);

[coeff_exec,score_exec]=pca(exec);


prepdata.FR=prep;
prepdata.ndist=ndist_prep;
prepdata.ndir=ndir_prep;
prepdata.npos=npos_prep;
prepdata.scores=score_prep;
prepdata.coeff=coeff_prep;


execdata.FR=exec;
execdata.ndist=ndist_exec;
execdata.ndir=ndir_exec;
execdata.npos=npos_exec;
execdata.ncycle=ncycle_exec;
execdata.scores=score_exec;
execdata.coeffs=coeff_exec;


end