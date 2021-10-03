function mink_sort = MinkSumEdgeSort2D(s1, s2, M1, M2)
% MinkSumEdgeSort2D Computes Minkowski sums between two convex bodies using
% edge sort method, works only for 2D
%
%  Inputs:
%    s1, s2   : Geometry class
%    M1, M2   : Linear transformation matrix
%
%  Outputs:
%    mink_sort: Minkowski sums points on the boundary
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021

dim = size(M1,1);
if dim ~= 2
    error('Only works for 2D!')
end

x1 = M1 * s1.GetPoints();
x2 = M2 * s2.GetPoints();

%% Compute Minkowski sum using edge sort algorithm
mink_sort = nan(dim,size(x1,2)+size(x2,2));

% the first vertex must be the lowest
x1 = reorder_polygon(x1);
x2 = reorder_polygon(-x2);

% we must ensure cyclic indexing
x1 = [x1, x1(:,1), x1(:,2)];
x2 = [x2, x2(:,1), x2(:,2)];

% main part
i = 1;
j = 1;
num = 1;

mink_sort(:,1) = x1(:,1) + x2(:,1);
while ~(i == size(x1,2)-1 && j == size(x2,2)-1)
    c = angle( x1(:,i+1)-x1(:,i), x2(:,j+1)-x2(:,j) );
    
    if i == size(x1,2)-1
        j = j+1;
        
    elseif j == size(x2,2)-1
        i = i+1;
        
    elseif c > 0
        i = i+1;
        
    elseif c <=0
        j = j+1;
        
    else
        i = i+1;
        j = j+1;
    end
    
    num = num+1;
    mink_sort(:,num) = x1(:,i) + x2(:,j);
end

end

%% Untility functions
function x_reordered = reorder_polygon(x)
% reorder_polygon Reorder polygon from minimum y-coordinates

[~, idx_min_y] = min(x(2,:));

if idx_min_y == 1
    x_reordered = x;
    return;
end

x_reordered = [x(:,idx_min_y:end), x(:,1:idx_min_y-1)];
end

function th = angle(p, q)
% angle Angle between two vectors p and q
th =  p(1) * q(2) - p(2) * q(1);
end