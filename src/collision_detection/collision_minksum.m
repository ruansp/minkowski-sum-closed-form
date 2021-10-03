function [status, dist, pt_cls, F_opt, m1_opt] = collision_minksum(s1, s2, opt)
% Solver for collision detection between two superquadrics
%
%  Inputs:
%    s1, s2: class of SuperQuadrics
%
%  Output:
%    status : status of contact
%               true -- in collision
%               false -- separated
%    dist   : Distance between the two objects
%    pt_cls : Witness points on s1, s2 and Minkowski sum boundary
%    F_opt  : Optimized cost function
%    m1_opt : Optimizer for gradient
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    MinkSumClosedForm, fsolve

R1 = quat2rotm(s1.q);
R2 = quat2rotm(s2.q);

% Initial guess
% Transform into local frame of s1
p0 = s2.tc - s1.tc;
s2_tc_in_s1 = R1' * p0;

% Optimization to solve min_{psi} |p0 - x_{1+2}(psi)|
% Initial condition
theta0 = [atan2( s2_tc_in_s1(3), norm(s2_tc_in_s1(1:2)) ),...
    atan2( s2_tc_in_s1(2), s2_tc_in_s1(1) )];

% Construct MinkSumClosedForm class object
minkObj = MinkSumClosedForm(s1, s2, R1, R2);

% Optimization
option = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt',...
    'Display', 'none',...
    'MaxIterations', 100);
[theta_opt, F_opt] = fsolve(@(theta) func(theta, minkObj, opt), theta0,...
    option);
m1_opt = minkObj.s1.GetGradientsFromSpherical(theta_opt);

% Distance and closest points
pt_cls.s1 = R1 * s1.GetPointsFromSpherical(theta_opt) + s1.tc;
pt_cls.mink = minkObj.GetMinkSumFromGradient(m1_opt) + s1.tc;
pt_cls.s2 = pt_cls.s1 + (s2.tc - pt_cls.mink);

dist = norm(pt_cls.s1 - pt_cls.s2);

% Collision status
status = true;
if strcmp(opt, 'ray')
    status = norm(pt_cls.s2 - s1.tc) < norm(pt_cls.mink - s1.tc);
elseif strcmp(opt, 'common-normal')
    status = (R1 * m1_opt)' * (pt_cls.s2 - pt_cls.s1) < 0;
end

if status; dist = -dist; end
    
end

function F = func(theta, minkObj, opt)
% Solver for collision detection between a point and the closed-form
% Minkowski sum boundary
%  Input:
%    theta: angular parameters for surface gradient
%
%  Output:
%    F    : equation for common normal constraint
%
% Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021

p0 = minkObj.s2.tc - minkObj.s1.tc;
m1 = minkObj.s1.GetGradientsFromSpherical(theta);
p_mink = minkObj.GetMinkSumFromGradient(m1);

if strcmp(opt, 'common-normal')
    % Common normal constraint
    m1_rot = quat2rotm(minkObj.s1.q) * m1;
    mink_to_p0 = p0 - p_mink;
    
    F = cross(m1_rot, mink_to_p0);
    
elseif strcmp(opt, 'ray')
    % Point-to-line distance
    F = cross(p0, p_mink);
end
end