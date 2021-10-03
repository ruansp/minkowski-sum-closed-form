% Demonstration on Minkowski sums of two convex bodies using definition
%  2D examples are shown. Algorithms include "convex hull", "edge sort" and
%  the propose "closed-form"
%
%  Author: 
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    SuperEllipse, MinkSumDefinition, MinkSumEdgeSort2D, MinkSumClosedForm

close all; clear; clc;
add_paths();

disp('***************************************************************************************')
disp('* Demonstration of Minkowski sums using convex hull and proposed closed-form solution *')
disp('***************************************************************************************')

mm = 50;
s1 = SuperEllipse([[10,8.6],0.36,0,0,0,0,mm]);
s2 = SuperEllipse([[9.4,6.8],1.28,0,0,0,0,mm]);

%% Transformations for the two SQs
angles = [pi/8,-pi/5];
shears = [1.32,1.6];
M1 = rot2(angles(1)) * [1, shears(1); 0, 1];
M2 = rot2(angles(2)) * [1, shears(2); 0, 1];

%% Minkowski sums
% Minkowski sums by definition
[mink_def_conv, mink_def, mink_k] = MinkSumDefinition(s1, s2, M1, M2);

% Minkowski sums by edge sort
mink_sort = MinkSumEdgeSort2D(s1, s2, M1, M2);

% Closed-form Minkowski sums
mink_closed_obj = MinkSumClosedForm(s1, s2, M1, M2);
m1 = s1.GetGradients();
mink_closed = mink_closed_obj.GetMinkSumFromGradient(m1);

%% Demo for definition
figure; hold on; axis equal; axis off;

% boundary of S1 and S2
m1 = s1.GetGradients();
x1 = M1 * s1.GetPointsFromGradient(m1);
patch(x1(1,:), x1(2,:)+25, 'r', 'LineWidth', 1.5);

m2 = s2.GetGradients();
x2 = M2 * s1.GetPointsFromGradient(m2);
patch(x2(1,:), x2(2,:)-25, 'b', 'LineWidth', 1.5);

% Points on the surfaces
x1 = M1 * s1.GetPointsFromGradient(m1);
plot(x1(1,:), x1(2,:)+25, 'k*', 'LineWidth', 1);

x2 = M2 * s1.GetPointsFromGradient(m2);
plot(x2(1,:), x2(2,:)-25, 'k*', 'LineWidth', 1);

% Points of S1+S2
plot(mink_def(1,:)+80, mink_def(2,:), 'k.', 'LineWidth', 1)

% Convex hull surface of S1+S2
plot(mink_def_conv(1,:)+180, mink_def_conv(2,:), ...
    'k*', 'LineWidth', 1)
patch(mink_def_conv(1,:)+180, mink_def_conv(2,:), ...
    'c', 'LineWidth', 1)

%% Demo for edge sort
figure; hold on; axis equal; axis off;

% boundary of S1 and S2
m1 = s1.GetGradients();
x1 = M1 * s1.GetPointsFromGradient(m1);
patch(x1(1,:), x1(2,:)+25, 'r', 'LineWidth', 1.5);

m2 = s2.GetGradients();
x2 = M2 * s1.GetPointsFromGradient(m2);
patch(x2(1,:), x2(2,:)-25, 'b', 'LineWidth', 1.5);

% Points on the surfaces
x1 = M1 * s1.GetPointsFromGradient(m1);
plot(x1(1,:), x1(2,:)+25, 'k*', 'LineWidth', 1);

x2 = M2 * s1.GetPointsFromGradient(m2);
plot(x2(1,:), x2(2,:)-25, 'k*', 'LineWidth', 1);

% Points of S1+S2
plot(mink_sort(1,:)+80, mink_sort(2,:), 'k*', 'LineWidth', 1)
patch(mink_sort(1,:)+80, mink_sort(2,:), 'c', 'LineWidth', 1)

%% Demo for closed-form
figure; hold on; axis equal; axis off;

% boundary of S1 and S2
m1 = s1.GetGradients();
x1 = M1 * s1.GetPointsFromGradient(m1);
patch(x1(1,:), x1(2,:)+25, 'r', 'LineWidth', 1.5);

m2 = s2.GetGradients();
x2 = M2 * s1.GetPointsFromGradient(m2);
patch(x2(1,:), x2(2,:)-25, 'b', 'LineWidth', 1.5);

% Points on the surface S1
x1 = M1 * s1.GetPointsFromGradient(m1);
plot(x1(1,:), x1(2,:)+25, 'k*', 'LineWidth', 1);

% Points of S1+S2
plot(mink_closed(1,:)+80, mink_closed(2,:), 'k*', 'LineWidth', 1)
patch(mink_closed(1,:)+80, mink_closed(2,:), 'c', 'LineWidth', 1)