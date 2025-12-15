function Mshared=sharedV_across_subspaces(PCrot_H,Data,idx_cond)
PCrot_H=PCrot_H./vecnorm(PCrot_H);
Nspaces=size(PCrot_H,3);
Mshared=nan(Nspaces,Nspaces);

%ndim=size(PCrot_H,2);

for i_space=1:Nspaces
    for j_space=i_space:Nspaces
        Mshared(i_space,j_space)=shared_variance_Co(Data(idx_cond==i_space,:),PCrot_H(:,:,i_space),PCrot_H(:,:,j_space));
        Mshared(j_space,i_space)=shared_variance_Co(Data(idx_cond==j_space,:),PCrot_H(:,:,j_space),PCrot_H(:,:,i_space));
        %         Mshared(i_space,j_space)=shared_variance(cov(Data(idx_cond==i_space,:)),PCrot_H(:,:,j_space),ndim);
        %         Mshared(j_space,i_space)=shared_variance(cov(Data(idx_cond==j_space,:)),PCrot_H(:,:,i_space),ndim);
    end
end
end

function V=shared_variance_Co(Data1,PC1,PC2)

%Just project and compute
score1=Data1*PC1;
OriVar=sum(sum(var(Data1)));
Var1=sum(var(score1))./OriVar;

score2=Data1*PC2;
Var2=sum(var(score2))./OriVar;

V=Var2/Var1;

end