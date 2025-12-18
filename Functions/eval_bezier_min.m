function error = eval_bezier_min(centre,P2,P1P3)
%EVAL_BEZIER_MIN Compute the fitting error between a quadratic Bézier curve
% and a set of reference points.
%
%   INPUTS:
%       centre - Reference points to be approximated by the Bézier curve.
%                This is an N-by-D matrix, where N is the number of points
%                (cycles) and D is the spatial dimension.
%
%       P2     - Second (middle) control point of the quadratic Bézier
%                curve, given as a 1-by-D vector.
%
%       P1P3   - First and third control points of the quadratic Bézier
%                curve, given as a 2-by-D matrix:
%                P1P3(1,:) is the first control point (P1),
%                P1P3(2,:) is the third control point (P3).
%
%   OUTPUT:
%       error  - Scalar value representing the total fitting error.
%                This is computed as the sum of the absolute differences
%                between each reference point and its closest point on
%                the Bézier curve.
%
%
% Andrea Colins Rodriguez
% 18/12/2025

Ncycles=size(centre,1);
s=0:0.05:1;

%% Evaluate Bezier with P2 guess
P(1,:)=P1P3(1,:);
P(2,:)=P2; 
P(3,:)=P1P3(2,:);

B = eval_bezier(P,s);

%% find minimal distance between curve and centre points (B_min)
B_min=nan(Ncycles,size(centre,2));

for i_cycle=1:Ncycles
    [~,idx_min]=min(pdist2(B,centre(i_cycle,:)));
    B_min(i_cycle,:)=B(idx_min,:);
end
error=sum(abs(B_min-centre),'all');
end