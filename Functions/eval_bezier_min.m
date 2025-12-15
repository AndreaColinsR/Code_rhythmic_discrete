function error=eval_bezier_min(center,P2,P1P3)
Ncycles=size(center,1);
t_in=0:0.05:1;
Nt=numel(t_in);
B=nan(Nt,size(center,2));

P1=P1P3(1,:);
P3=P1P3(2,:);

for ti=1:Nt
    t=t_in(ti);
    B(ti,:)=(1-t)*(P1+t*P2)+t*((1-t)*P2+t*P3);
end
B_min=nan(Ncycles,size(center,2));

%% find points with minimal distance (B_min)
for i_cycle=1:Ncycles
    [~,idx_min]=min(pdist2(B,center(i_cycle,:)));
    B_min(i_cycle,:)=B(idx_min,:);
end
error=sum(abs(B_min-center),'all');
end