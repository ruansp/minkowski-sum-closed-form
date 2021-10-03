function [num_contacts, dist] = avg_dist_implicit(s1, s2, M1, M2, mink)
% avg_dist_implicit Averaged distance between point set and implicit
% function of a superquadric model
%
%  Inputs:
%    s1, s2       : Geometry classes
%    M1, M2       : Linear tranformation matrices
%    mink         : A set of points on Minkowski sums boundary
%
%  Output:
%    num_contacts : Number of contact points for the point set of 
%                   Minkowski sums boundary
%    dist         : Averaged distance of point set to implicit function
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    implicit_function

is_intersect = zeros(size(mink,2),1);
min_dist = zeros(size(mink,2),1);

for i = 1:size(mink,2)
    % distances of points, from implicit function
    val = implicit_function(s1,s2,M1,M2,mink(:,i));
    
    % indicate intersections
    is_intersect(i) = sum(val < 0) > 0;
    
    % compute minimum value as distance measure
    min_dist(i) = min(val);
end

dist = sum(min_dist)/size(mink,2);
num_contacts = sum(is_intersect)/size(mink,2);

% Remove outliers
if dist > 1
    dist = 0;
    num_contacts = 0;
end