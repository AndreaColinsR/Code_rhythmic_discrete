function [EigsA,r2,r1,A] = LDS(scores)
% LDS Fits a Linear Dynamical System (LDS) to a low-dimensional trajectory segment.
% The dynamics matrix A is estimated using a least-squares approximation
%   of the temporal derivative computed via a central difference scheme.
%   If the number of time points is insufficient (T ≤ D²), the function
%   returns NaN outputs.
%
%   [EigsA,r2,r1,A] = LDS(scores)
%
%   The function estimates a linear dynamical system of the form
%       dy/dt = A y
%   from a low-dimensional time series and evaluates how well the fitted
%   dynamics reconstruct the original trajectory.
%
%   INPUTS
%   ------
%
%   scores : [T x D] matrix
%       Low-dimensional trajectory (e.g. PCA scores), where T is the number
%       of time points and D is the dimensionality of the state space.
%
%   OUTPUTS
%   -------
%
%   EigsA : [D x 1] vector
%       Eigenvalues of the fitted dynamics matrix A, characterising the
%       stability and oscillatory properties of the system.
%
%   r2 : scalar
%      Squared (Pearson) correlation between the
%       reconstructed trajectory and the original data.
%
%   r1 : [T x D] matrix
%       Reconstructed trajectory obtained by integrating the fitted linear
%       dynamical system forwards and backwards in time from the midpoint.
%
%   A : [D x D] matrix
%       Estimated linear dynamics matrix defining the system dy/dt = A y.
%
% Andrea Colins
% 22/12/2025

debugging=0;

npoints=size(scores,1);
ndim=size(scores,2);
EigsA=nan(ndim,1);
r2=nan;
r1=nan;

% check that scores has enough data points to fit the model
if npoints<=ndim^2
    disp('The trajectory has very few timebins (Timebins<=Dim^2)')
    return
end


% approximate derivative as the central difference (one bin forward and one bin backwards)
  dy = diff(scores);
  A  = (((scores(2:end,:)+scores(1:end-1,:))/2)\dy)';

EigsA=eig(A);

%% evaluate the fitted curve
r1=scores';
middlepoint=round(npoints/2);

%% Eq: dy/dt=A*y
%forward
for i=middlepoint:npoints
    
    r1(:,i)= A*r1(:,i-1)+r1(:,i-1);
end

%backward
for i=middlepoint-1:-1:1
    
    r1(:,i)= -A*r1(:,i+1)+r1(:,i+1);
end

r1=r1';

%% Evaluate performance
r2=corr(r1(:),scores(:)).^2;


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