% Demonstration on two theorems of computing closed-form Minkowski sums of
% two convex bodies
%
%  Author: 
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    SuperEllipse, MinkSumClosedForm

close all; clear; clc;
add_paths();

disp('*************************************************************')
disp('* Demonstration of the proposed closed-form Minkowski sums  *')
disp('*************************************************************')

%% Original shapes of two SQs
mm = 100;
s1 = SuperEllipse([[3.3,2.6],0.73,0,0,0,0,mm]);
s2 = SuperEllipse([[2.4,2],1.28,0,0,0,0,mm]);

angles = [pi/6.3, -pi/3];
shears = [0.276, 0.85];
M1 = rot2(angles(1)) * [1, shears(1); 0, 1];
M2 = rot2(angles(2)) * [1, shears(2); 0, 1];

x1 = M1 * s1.GetPointsCanonical();
x2 = M2 * s2.GetPointsCanonical();

%% Closed-form Minkowski sums
idx = 25;
mink_closed_obj = MinkSumClosedForm(s1, s2, M1, M2);
m1 = s1.GetGradients();
mink_closed = mink_closed_obj.GetMinkSumFromGradient(m1);
x_mink = mink_closed(:,idx);
x2 = x2 + x_mink;

%% Contact parameters
m1 = M1' \ s1.GetGradients();
m1 = m1(:,idx);
n1 = m1/norm(m1);

phi_m1 = mink_closed_obj.Phi(m1);
m2 = phi_m1/norm(m1) * m1;

%% First theorem (Theorem 4.1 in the article)
disp('** Normal parameterization')

figure; hold on; axis equal; axis off;
% Original surfaces
patch(x1(1,:), x1(2,:), 'g', 'FaceAlpha', 0.3)
patch(x2(1,:), x2(2,:), 'b', 'FaceAlpha', 0.3)

% Body centers
plot(0, 0, '*k', 'LineWidth', 2)
plot(x_mink(1), x_mink(2), '*k', 'LineWidth', 2)

% Contact point and normals
plot(x1(1,idx), x1(2,idx), '*k', 'LineWidth', 2)
arrow3(x1(:,idx)', x1(:,idx)'+n1', '-k1.5', 0.8, 1.5)
arrow3(x1(:,idx)', x1(:,idx)'-n1', '-k1.5', 0.8, 1.5)

arrow3([0 0], x1(:,idx)', '--k1', 0.6, 1.5)
arrow3(x_mink', x1(:,idx)', '--k1', 0.6, 1.5)
arrow3([0 0], x_mink', '--k1', 0.6, 1.5)

% Part of Mink sum boundary
num_neighbors = 20;
plot(mink_closed(1,idx-num_neighbors:idx+num_neighbors),...
    mink_closed(2,idx-num_neighbors:idx+num_neighbors),...
    'r', 'LineWidth', 3)

%% Second theorem (Theorem 4.3 in the article)
disp('** Gradient (un-normalized) parameterization')

figure; hold on; axis equal; axis off;
% Original surfaces
patch(x1(1,:), x1(2,:), 'g', 'FaceAlpha', 0.3)
patch(x2(1,:), x2(2,:), 'b', 'FaceAlpha', 0.3)

% Body centers
plot(0, 0, '*k', 'LineWidth', 2)
plot(x_mink(1), x_mink(2), '*k', 'LineWidth', 2)

% Contact point and normals
plot(x1(1,idx), x1(2,idx), '*k', 'LineWidth', 2)
arrow3(x1(:,idx)', x1(:,idx)'+m1', '-k1.5', 0.8, 1.5)
arrow3(x1(:,idx)', x1(:,idx)'-m2', '-k1.5', 0.8, 1.5)

arrow3([0 0], x1(:,idx)', '--k1', 0.6, 1.5)
arrow3(x_mink', x1(:,idx)', '--k1', 0.6, 1.5)
arrow3([0 0], x_mink', '--k1', 0.6, 1.5)

% Part of Mink sum boundary
num_neighbors = 20;
plot(mink_closed(1,idx-num_neighbors:idx+num_neighbors),...
    mink_closed(2,idx-num_neighbors:idx+num_neighbors),...
    'r', 'LineWidth', 3)