function [execdata,prepdata,endexecdata,execdata_all]=get_prep_exec_after_FR(FR,idx_pos,idx_dir,idx_dist,idx_Ncycle)


Ndir=max(idx_dir);
Npos=max(idx_pos);
Ncycle=max(idx_Ncycle);
Ndist=unique(idx_dist);

prep=[];
exec=[];
end_exec=[];
exec_all=[];


ndist_prep=[];
ndir_prep=[];
npos_prep=[];



ncycle_exec=[];
ndist_exec=[];
ndir_exec=[];
npos_exec=[];
npos_endexec=[];

ncycle_exec_all=[];
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
            exec_all=[exec_all;FR(cond(idx_onset+1:idx_mov_end),:)];
            ndist_exec_all=[ndist_exec_all;ones(Nexec,1)*i_dist];
            ndir_exec_all=[ndir_exec_all;ones(Nexec,1)*i_dir];
            npos_exec_all=[npos_exec_all;ones(Nexec,1)*i_pos];
            
            %% for each condition select exec of first half cycle for dpca
            % start of second cycle or end of sequence
            if Ndist(i_dist)<1
                % for half a cycle, take the entire movement
                idx_offset=idx_mov_end;
                idx_onset3=idx_onset;
                    
                % or
                % for half a cycle, take only the half that is
                % accelerating/desaccelerating. 

                %idx_offset=round((idx_mov_end+idx_onset)/2);
                %idx_onset3=round((idx_onset+idx_offset)./2);
               
            elseif Ndist(i_dist)==1
                 % for one cycle, take the first half cycle
                firstcycledur=find(isnan(idx_Ncycle(cond(idx_onset+1:end))),1,"first");
                idx_offset=idx_onset+round(firstcycledur/2);
                
                % for one cycle, take the second half cycle
                idx_onset3=idx_offset;

            else
                % for more than one cycle, take the first half cycle
                firstcycledur=find(idx_Ncycle(cond(idx_onset+1:end))>1,1,"first");
                idx_offset=idx_onset+round(firstcycledur/2);

                % for one cycle, take the second half cycle
                startcycledur=find(idx_Ncycle(cond(idx_onset+1:end))<Ndist(i_dist),1,"last");
                lastcycledur=idx_mov_end-idx_onset+1-startcycledur;
                idx_onset3=idx_onset+startcycledur+round(lastcycledur./2);
            end
            
            % to try same segment for all
             %idx_offset=idx_onset2+300;
             %idx_onset3=idx_mov_end-300;

            FRfirstHalf=interp1(linspace(0,1,idx_offset-idx_onset2+1)',FR(cond(idx_onset2:idx_offset),:),linspace(0,1,tmp)');
            % last cycle
            %FRlastHalf=interp1(linspace(0,1,idx_mov_end-idx_onset3+1)',FR(cond(idx_onset3:idx_mov_end),:),linspace(0,1,tmp)');
            
            % After mov
            FRlastHalf=interp1(linspace(0,1,numel(cond(idx_mov_end:end)))',FR(cond(idx_mov_end:end),:),linspace(0,1,tmp)');
            
          

            exec=[exec;FRfirstHalf];
            end_exec=[end_exec;FRlastHalf];

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

%baseline=mean(baseline(:,selected_neurons),1);

% soft normalise
norm_factor=range(FR)+5;
prep=prep(:,selected_neurons)./norm_factor;
exec=exec(:,selected_neurons)./norm_factor;
end_exec=end_exec(:,selected_neurons)./norm_factor;
exec_all=exec_all(:,selected_neurons)./norm_factor;

[coeff_prep,score_prep]=pca(prep);

[coeff_exec,score_exec]=pca(exec);

[coeff_endexec,score_endexec]=pca(end_exec);

[coeff_exec_all,score_exec_all]=pca(exec_all);

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


endexecdata.FR=end_exec;
endexecdata.ndist=ndist_exec;
endexecdata.ndir=ndir_exec;
endexecdata.npos=npos_endexec;
endexecdata.ncycle=ncycle_exec;
endexecdata.scores=score_endexec;
endexecdata.coeffs=coeff_endexec;

execdata_all.FR=exec_all;
execdata_all.ndist=ndist_exec_all;
execdata_all.ndir=ndir_exec_all;
execdata_all.npos=npos_exec_all;
execdata_all.scores=score_exec_all;
execdata_all.coeffs=coeff_exec_all;
end