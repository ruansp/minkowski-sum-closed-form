% Demo script for exact closed-form Minkowski sums between two 2D
% superellipses
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    SuperEllipse, MinkSumClosedForm

close all; clear; clc;
add_paths();

disp('*********************************************************************')
disp('* Demonstration of closed-form Minkowski sums for two superellipses *')
disp('*********************************************************************')

s1 = SuperEllipse([1+5*rand(1,2), 0.1+1.8*rand, -0.5+rand, 0, 0, 0, 100]);
s2 = SuperEllipse([1+5*rand(1,2), 0.1+1.8*rand, -0.5+rand, 0, 0, 0, 100]);

% Transformations for the two SQs
angles = 2*pi*rand(1,2);
shears = 2*(2*rand(1,2)-1);

M1 = rot2(angles(1)) * [1, shears(1); 0, 1];
M2 = rot2(angles(2)) * [1, shears(2); 0, 1];

% Compute Minkowski sums
mink_obj = MinkSumClosedForm(s1, s2, M1, M2);
m1 = s1.GetGradientsCanonical();
mink = mink_obj.GetMinkSumFromGradient(m1);

%% Plots
figure; hold on; axis equal; axis off;
% Plot Minkowski sum
plot(mink(1,:), mink(2,:), 'r*', 'LineWidth', 3);

% Plot s1
x1 = M1 * s1.GetPointsFromGradient(m1);
patch(x1(1,:), x1(2,:), 'k', 'LineWidth', 0.1, 'FaceAlpha', 0.3);

% Plot s2
for i = 1:size(mink,2)
    s2_touch = M2 * s2.GetPointsCanonical() + mink(:,i);
    plot(s2_touch(1,:), s2_touch(2,:), 'b', 'LineWidth', 0.1);
end

%% Demo on sampled s2
figure; hold on; axis equal; axis off;
% Plot Minkowski sum
plot(mink(1,:), mink(2,:), 'r*', 'LineWidth', 3);

% Plot s1
x1 = M1 * s1.GetPointsFromGradient(m1);
patch(x1(1,:), x1(2,:), 'k', 'LineWidth', 0.1, 'Facealpha', 0.3);

% Plot s2
k = 3;
s2_touch = M2 * s2.GetPointsCanonical() + mink(:,i);
plot(s2_touch(1,:), s2_touch(2,:), 'b.', 'LineWidth', 0.1);