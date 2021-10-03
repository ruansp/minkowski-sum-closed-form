% Main file for application on collision detection between two bodies and
% compared to common normal method. Demonstration of ideas
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    SuperQuadrics, Ellipsoid, collision_minksum, collision_minksum, 
%    collision_parametric

close all; clear; clc;

add_paths();
mm = [100, 100];

s1 = SuperQuadrics({[6.4, 4.35, 3.67], [0.8, 1.6], [0,0],...
    [-2.4, 2.5, -1.3]', [0.2, 0.13, 0.2, 0.89], mm});
s2 = SuperQuadrics({[3.76, 8.33, 2.64], [1.434, 0.56], [0,0],...
    [-12.4, 23.54, 3.46]', [-0.3, 0.23, -0.2, 0.92], mm});

%% Minkowski sums
minkObj = MinkSumClosedForm(s1, s2, quat2rotm(s1.q), quat2rotm(s2.q));
m1 = s1.GetGradientsCanonical();

mink = minkObj.GetMinkSumFromGradient(m1) + s1.tc;
mink_x = reshape(mink(1,:), mm);
mink_y = reshape(mink(2,:), mm);
mink_z = reshape(mink(3,:), mm);

% Minkowski-based (ray) status
[status_mink_ray, dist_mink_ray, pt_mink_ray, f_mink_ray, m_mink_ray] =...
    collision_minksum(s1, s2, 'ray');

% Minkowski-based (common normal) status
[status_mink_normal, dist_mink_normal, pt_mink_normal, f_mink_normal, m_mink_normal] =...
    collision_minksum(s1, s2, 'common-normal');

% Param surface common normal status
[status_param, dist_param, pt_param, f_param, theta_param] =...
    collision_parametric(s1, s2);

%% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Common normal method %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; axis equal; axis off;
lightangle(gca,90,30);
lighting gouraud;
view([363,100,10])

s1.PlotShape('g', 1);
s2.PlotShape('g', 1);

% Witness points
plot3(pt_param.s1(1), pt_param.s1(2), pt_param.s1(3),...
    'k*', 'LineWidth', 2)
plot3(pt_param.s2(1), pt_param.s2(2), pt_param.s2(3),...
    'k*', 'LineWidth', 2)
plot3([pt_param.s1(1) pt_param.s2(1)],...
    [pt_param.s1(2) pt_param.s2(2)],...
    [pt_param.s1(3) pt_param.s2(3)], 'k--', 'LineWidth', 1)

% Solved normal vectors
m1_param_opt = s1.GetGradientsFromSpherical(theta_param(1:2));
m2_param_opt = s2.GetGradientsFromSpherical(theta_param(3:4));
arrow3([pt_param.s1, pt_param.s2]',...
    [(pt_param.s1 +...
    3*quat2rotm(s1.q)*(m1_param_opt/norm(m1_param_opt))),...
    (pt_param.s2 +...
    3*quat2rotm(s2.q)*(m2_param_opt/norm(m2_param_opt)))]',...
    'k-1.5', 0.6, 1)

%%%%%%%%%%%%%%%%%%%%%%%
%%% Minkowski (Ray) %%%
%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; axis equal; axis off;
lightangle(gca,90,30);
lighting gouraud;
view([363,100,10])

s1.PlotShape('g', 0.5);
s2.PlotShape('g', 0.5);

% Center points
plot3(s1.tc(1), s1.tc(2), s1.tc(3), 'b*', 'LineWidth', 2)
plot3(s2.tc(1), s2.tc(2), s2.tc(3), 'b*', 'LineWidth', 2)
plot3([s1.tc(1) s2.tc(1)], [s1.tc(2) s2.tc(2)], [s1.tc(3) s2.tc(3)],...
    'b-', 'LineWidth', 2)

% Solved point on Minkowski sums boundary
plot3(pt_mink_ray.mink(1), pt_mink_ray.mink(2),...
    pt_mink_ray.mink(3), 'r*', 'LineWidth', 2)
surf(mink_x, mink_y, mink_z,...
    'FaceColor', 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Minkowski (Common normal) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; axis equal; axis off;
lightangle(gca,90,30);
lighting gouraud;
view([363,100,10])

s1.PlotShape('g', 0.5);
s2.PlotShape('g', 0.5);

% Center points
plot3(s2.tc(1), s2.tc(2), s2.tc(3), 'b*', 'LineWidth', 2)

% Solved point on Minkowski sums boundary
plot3(pt_mink_normal.mink(1), pt_mink_normal.mink(2),...
    pt_mink_normal.mink(3), 'r*', 'LineWidth', 2)
surf(mink_x, mink_y, mink_z,...
    'FaceColor', 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
plot3([pt_mink_normal.mink(1) s2.tc(1)],...
    [pt_mink_normal.mink(2) s2.tc(2)],...
    [pt_mink_normal.mink(3) s2.tc(3)],...
    'b--', 'LineWidth', 1)
arrow3(pt_mink_normal.mink',...
    ( pt_mink_normal.mink +...
    5*quat2rotm(s1.q)*(m_mink_normal/norm(m_mink_normal)) )',...
    '-k1.5', 0.6, 1.5)

% Witness points
plot3(pt_mink_normal.s1(1), pt_mink_normal.s1(2),...
    pt_mink_normal.s1(3), 'k*', 'LineWidth', 2)
plot3(pt_mink_normal.s2(1), pt_mink_normal.s2(2),...
    pt_mink_normal.s2(3), 'k*', 'LineWidth', 2)
plot3([pt_mink_normal.s1(1) pt_mink_normal.s2(1)],...
    [pt_mink_normal.s1(2) pt_mink_normal.s2(2)],...
    [pt_mink_normal.s1(3) pt_mink_normal.s2(3)], 'k--', 'LineWidth', 1)
