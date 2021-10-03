function [mink_def_conv, mink_def, k] = MinkSumDefinition(s1, s2, M1, M2)
% MinkSumDefinition Computes Minkowski sums between two convex bodies using
% original definition and convex hull
%
%  Inputs:
%    s1, s2       : Geometry class
%    M1, M2       : Linear transformation matrix
%
%  Outputs:
%    mink_def_conv: The convex hull of the Minkowski sums
%    mink_def     : Direct additions using Minkowski sums definition
%    k            : Convex hull in terms of a vector of point indices 
%                   arranged in a counter-clockwise cycle around the hull.
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    convhull

dim = size(M1,1);

x1 = M1 * s1.GetPoints();
x2 = M2 * s2.GetPoints();

%% Compute Minkowski sum point-by-point from definition
mink_def = zeros(dim,size(x1,2)*size(x2,2));
% Element-wise addition
for i = 1:dim
    temp = x1(i,:) + x2(i,:)';
    mink_def(i,:) = temp(:);
end

% Compute convex hull
if dim == 2
    k = convhull(mink_def(1,:), mink_def(2,:));
elseif dim ==3
    k = convhull(mink_def(1,:), mink_def(2,:), mink_def(3,:));
end
mink_def_conv = mink_def(:,k);
end