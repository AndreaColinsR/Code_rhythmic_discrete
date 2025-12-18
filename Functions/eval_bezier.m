function B = eval_bezier(P,s)
%EVAL_BEZIER2 Evaluate a quadratic Bézier curve using the
%       de Casteljau algorithm. 
%
%
%   INPUTS:
%       P     - Control points of the Bézier curve, given as a 3-by-D
%               matrix, where each row represents one control point and
%               D is the spatial dimension (e.g., D = 2 for 2D, D = 3 for 3D).
%               P(1,:) is the first control point,
%               P(2,:) is the second control point,
%               P(3,:) is the third control point.
%
%       s  - Vector of parameter values at which to evaluate the curve.
%               Each value should typically satisfy 0 <= s <= 1.
%
%   OUTPUT:
%       B     - Evaluated points on the Bézier curve. This is an
%               length(t_in)-by-D matrix, where each row corresponds to
%               the Bézier curve point evaluated at the corresponding
%               parameter value in t_in.
%
% Andrea Colins Rodriguez
% 18/12/2025

Nt = numel(s);
B = nan(Nt,size(P,2));

for ti = 1:Nt
    t = s(ti);
    B(ti,:) = (1-t)*(P(1,:)+t*P(2,:))+t*((1-t)*P(2,:)+t*P(3,:));
end

end