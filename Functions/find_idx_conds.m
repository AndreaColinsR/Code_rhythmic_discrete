function [idx_dir,idx_pos,idx_Ncycle,idx_dist,prep,exec,startexec,endexec]=find_idx_conds(trials_idx,idx_conds,idx_current_cycle)

Nconds=max(trials_idx);
ntimepoints=size(trials_idx,1);
idx_dir=zeros(ntimepoints,1);
idx_dist=zeros(ntimepoints,1);
idx_pos=zeros(ntimepoints,1);
idx_Ncycle=nan(ntimepoints,1);
exec=ones(ntimepoints,1);
prep=ones(ntimepoints,1);
endexec=ones(ntimepoints,1);
startexec=zeros(ntimepoints,1);
%NumberCyles=unique(idx_conds(:,1));



for i_cond=1:Nconds
    this_cond=find(trials_idx==i_cond);

    idx_Ncycle(this_cond)=idx_current_cycle(i_cond,1:numel(this_cond));
    idx_pos(this_cond)=idx_conds(i_cond,2);
    idx_dir(this_cond)=idx_conds(i_cond,3);

    exec(this_cond(1:100))=0;
    exec(this_cond(end-40:end))=0;
    % prep idx
    prep(this_cond(101:end))=0;


    idx_mov_end=find(~isnan(idx_current_cycle(i_cond,:)),1,'last');

    if idx_conds(i_cond,1)<1

        end_first_half=idx_mov_end;

    elseif idx_conds(i_cond,1)==1
        firstcycledur=idx_mov_end-100;
        end_first_half=100+round(firstcycledur/2);
    else
        firstcycledur=find(idx_current_cycle(i_cond,100:end)>1,1,"first");
        end_first_half=100+round(firstcycledur/2);

    end
    startexec(this_cond(100:end_first_half))=1;

    % last half cycle
    endexec(this_cond(1:end-100))=0;
    endexec(this_cond(end-40:end))=0;

    %idx_dist(this_cond)=find(idx_conds(i_cond,1)==NumberCyles);
    idx_dist(this_cond)=idx_conds(i_cond,1);
end

end