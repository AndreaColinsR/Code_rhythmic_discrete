function Init_cond_t=initial_cond_helix(states,idx_dir,idx_pos,idx_dist,idx_Ncycle,exec,do_plot)

exec=exec>0;
states=states(exec>0,:);
idx_dir=idx_dir(exec,:);
idx_pos=idx_pos(exec,:);
idx_dist=idx_dist(exec,:);
idx_Ncycle=idx_Ncycle(exec,:);

Ndist=unique(idx_dist,'stable');
[~,scores]=pca(states);


dt=1;
if size(states,1)>10000
    dt=10;
end

ndims=6;
scores=scores(:,1:ndims);
center=nan(7,ndims);

Init_cond_t=nan(2*2,numel(Ndist));

[dist_sorted,~]=sort(Ndist,'ascend');

t=0:0.01:1;

counter=1;

if do_plot
    colour_cycle=colormap(parula(7));
    colour_dist=plasma(numel(Ndist));
end

for i_dir=1:2
    for i_pos=1:2


        for i_cycle=1:7

            this_cond=find(idx_dir==i_dir & idx_pos==i_pos & idx_dist==7 & idx_Ncycle==i_cycle);

            center(i_cycle,:)=mean(scores(this_cond,1:ndims));
            if do_plot && i_dir==2 && i_pos==2
                plot3(scores(this_cond,1),scores(this_cond,2),scores(this_cond,3),'Color',colour_cycle(i_cycle,:))
                hold on
                plot3(center(i_cycle,1),center(i_cycle,2),center(i_cycle,3),'o','MarkerFaceColor',[0 0 0]+0.1*i_cycle,'MarkerEdgeColor',[0 0 0]+0.1*i_cycle,'MarkerSize',8)
            end
        end

        %% fit to bezier curve
        P(1,:)=center(1,:);
        P(2,:)=center(6,:);
        P(3,:)=center(end,:);

        P1P3=[P(1,:);P(3,:)];
        P2=P(2,:);
        %% order of parameters fcn = @(coefficients,problemparameters,x,y) expression
        ftmp = @(P2) eval_bezier_min(center,P2,P1P3);
        x = fminsearch(ftmp,P2);
        P(2,:)=x;
        B=eval_bezier2(P,t);

        if do_plot && i_dir==2 && i_pos==2
            plot3(B(:,1),B(:,2),B(:,3),'Color',[0.5 0.5 0.5])
        end


        InitC=nan(numel(Ndist),ndims);
        idx_min=nan(numel(Ndist),1);

        this_cond=idx_dir==i_dir & idx_pos==i_pos & idx_dist==7;

        %find the closest point 
        for i_dist=1:numel(Ndist)
            this_cond2=find(idx_dir==i_dir & idx_pos==i_pos & idx_dist==dist_sorted(i_dist));
            InitC(i_dist,:)=mean(scores(this_cond2(1:10*dt),:),1);
            [~,idx_min(i_dist)]=min(pdist2(scores(this_cond,:),InitC(i_dist,:)));

            if do_plot && i_dir==2 && i_pos==2
                plot3(scores(this_cond2,1),scores(this_cond2,2),scores(this_cond2,3),'Color',colour_dist(i_dist,:))
                hold on
                plot3(scores(this_cond2(1),1),scores(this_cond2(1),2),scores(this_cond2(1),3),'.','Color',colour_dist(i_dist,:),'MarkerSize',14)

            end

        end


        Init_cond_t(counter,:)=idx_min./sum(this_cond);%t(idx_min);

        counter=counter+1;
    end
end

%Init_cond_t=Init_cond_t-Init_cond_t(:,end);
end

function B=eval_bezier2(P,t_in)
Nt=numel(t_in);
B=nan(Nt,size(P,2));
for ti=1:Nt
    t=t_in(ti);
    B(ti,:)=(1-t)*(P(1,:)+t*P(2,:))+t*((1-t)*P(2,:)+t*P(3,:));
end
end

