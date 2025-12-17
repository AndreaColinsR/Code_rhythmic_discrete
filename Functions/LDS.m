function [EigsA,r2,r1,A]=LDS(scores)
%% LDS fits a Linear Dynamical System to a trajectory scores
% scores:
debugging=0;

npoints=size(scores,1);
ndim=size(scores,2);
EigsA=nan(ndim,1);
r2=nan;
r1=nan;

% check that scores has enough data points to fit the model
if npoints<=ndim^2
    return
end


% fit A*y=dy/dt

% 1) from Mark's code
% dy=diff(scores);
% A=(scores(1:end-1,:)\dy)';
% A=(scores(2:end,:)\dy)';

% 2) take the middle point between bins
  dy=diff(scores);
  A=(((scores(2:end,:)+scores(1:end-1,:))/2)\dy)';

% 3) take longer derivative
% dy=(scores(3:end,:)-scores(1:end-2,:))/2;
% A=(scores(2:end-1,:)\dy)';

EigsA=eig(A);
% EigsA2=eigs(A2);
% nim(nseg)=sum(abs(imag(EigsA))>0);

%% evaluate the fitted curve
r1=scores';
middlepoint=round(npoints/2);

%% Eq: dy/dt=A*y
%forward
for i=middlepoint:npoints
    
    % y(t)=A*y(t-1)*(t-t-1)+y(t-1)
    r1(:,i)= A*r1(:,i-1)+r1(:,i-1);
end

%backward
for i=middlepoint-1:-1:1
    
    r1(:,i)= -A*r1(:,i+1)+r1(:,i+1);
end

r1=r1';

r2=corr(r1(:),scores(:)).^2;


%% Evaluate performance
% debugging

if debugging
    subplot(2,2,3)
    plot3(scores(:,1),scores(:,2),scores(:,3))
    hold on
    plot3(r1(:,1),r1(:,2),r1(:,3))
    legend('Original','fitting')
    hold off
    
    subplot(2,2,4)
    hold on
    plot(scores(:))
    plot(r1(:))
    hold off
end

end