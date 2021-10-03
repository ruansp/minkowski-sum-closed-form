function [shape, cost] = sq_fitting(var_init, pts)
% sq_fitting Superquadric model fitting for point cloud data
%
%  Inputs:
%    var_init: inital values of variables
%              var_init(1:3) -- semi-axis lengths
%              var_init(4:5) -- exponents
%              var_init(6:8) -- position of the center
%              var_init(1:3) -- orientation
%    pts     : input point cloud data points to be fitted
%
%  Outputs:
%    shape   : Information of fitted superquadric shape
%              shape.a -- semi-axis lengths
%              shape.eps -- exponents
%              shape.tc -- position of the center
%              shape.q-- orientation
%    cost    : Optimized cost value of the fitting optimization
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    eig, lsqnonlin

if isempty(var_init)
    % Initial condition using PCA
    tc = mean(pts, 2);
    
    [U,S,~] = eig((pts-tc)*(pts-tc)');
    a = diag(S).^(0.5)/20;
    
    if det(U) < 0
        U(:,2) = -U(:,2);
    end
    q = rotm2quat(U);
    
    var0 = [a', [1,1], tc', q];
    
else
    var0 = var_init;
end

lb = [0,0,0,0.1,0.1,-inf,-inf,-inf,-1,-1,-1,-1];
ub = [inf,inf,inf,2,2,inf,inf,inf,1,1,1,1];

options = optimoptions('lsqnonlin',...
    'Algorithm', 'levenberg-marquardt', 'display', 'off');
[var_opt, cost] = lsqnonlin(@(var) sq_model(var, pts),...
    var0, lb, ub, options);

% Retrieve results
shape.a = var_opt(1:3);
shape.eps = var_opt(4:5);
shape.pose.tc = var_opt(6:8)';
shape.pose.q = var_opt(9:12);
end

%% Objective function to be minimized
function F = sq_model(var, pts)
% parameters
a1 = var(1);
a2 = var(2);
a3 = var(3);
eps1 = var(4);
eps2 = var(5);
tc = var(6:8)';
q = var(9:12) ./ norm(var(9:12));

% Transformation to the body frame of the superquadric
R = quat2rotm(q);
pts = R\(pts - tc);

% Cost function
f = eps_fun( eps_fun( (pts(1,:) ./ a1).^2, 1/eps2 ) +...
    eps_fun( (pts(2,:) ./ a2).^2, 1/eps2 ), eps2/eps1 ) +...
    eps_fun( (pts(3,:) ./ a3).^2, 1/eps1 );
F = sqrt(a1*a2*a3) * abs(f .^ eps1 -1);
end