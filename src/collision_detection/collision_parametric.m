function [status, dist, pt_cls, F_opt, theta_opt] = collision_parametric(s1, s2)
% collision_parametric answer proximity queries between two superquadrics,
% using parametric functions and common normal concept
%
%  Inputs:
%    s1, s2      : SuperQuadrics class
%
%  Output:
%    status    : status of contact
%                  true -- in collision
%                  false -- separated
%    dist      : Distance between the two objects
%    pt_cls    : Witness points on s1 and s2
%    F_opt     : Optimized cost function
%    theta_opt : Optimizer for angular parameters
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    fsolve

R1 = quat2rotm(s1.q);
R2 = quat2rotm(s2.q);

% Initial guess
% Transform into local frame of s1
p0 = s2.tc - s1.tc;
s2_tc_in_s1 = R1' * p0;
s1_tc_in_s2 = R2' * (-p0);

% Optimization to solve min_{psi} |p0 - x_{1+2}(psi)|
% Initial condition
theta0 = [atan2( s2_tc_in_s1(3), norm(s2_tc_in_s1(1:2)) ),...
    atan2( s2_tc_in_s1(2), s2_tc_in_s1(1) ),...
    atan2( s1_tc_in_s2(3), norm(s1_tc_in_s2(1:2)) ),...
    atan2( s1_tc_in_s2(2), s1_tc_in_s2(1) )];

% Optimization
option = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt',...
    'Display', 'none',...
    'MaxIterations', 100);
[theta_opt, F_opt] = fsolve(@(th) func(th, s1, s2), theta0, option);

% Distance and closest points
pt_cls.s1 = R1 * s1.GetPointsFromSpherical(theta_opt(1:2)) + s1.tc;
pt_cls.s2 = R2 * s2.GetPointsFromSpherical(theta_opt(3:4)) + s2.tc;
dist = norm(pt_cls.s2 - pt_cls.s1);

% Collision status
m1 = R1 * s1.GetGradientsFromSpherical(theta_opt(1:2));

if m1' * (pt_cls.s2 - pt_cls.s1) > 0
    status = 0;
else
    status = 1;
end
end

%% Objective function
function F = func(th, s1, s2)
R1 = quat2rotm(s1.q);
R2 = quat2rotm(s2.q);

m1 = s1.GetGradientsFromSpherical(th(1:2));
n1 = R1 * m1/norm(m1);
m2 = s2.GetGradientsFromSpherical(th(3:4));
n2 = R2 * m2/norm(m2);

x1 = R1 * s1.GetPointsFromSpherical(th(1:2)) + s1.tc;
x2 = R2 * s2.GetPointsFromSpherical(th(3:4)) + s2.tc;
x1_to_x2 = x2-x1;

% n1 x n2 + n1 x r12
F = [cross(n1, n2);
     cross(n1, x1_to_x2)];
end