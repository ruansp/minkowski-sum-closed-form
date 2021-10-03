% Demo script for exact closed-form Minkowski sums between two 3D
% superquadrics
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    SuperQuadrics, MinkSumClosedForm

close all; clear; clc;
add_paths();

disp('*********************************************************************')
disp('* Demonstration of closed-form Minkowski sums for two superquadrics *')
disp('*********************************************************************')

mm = [100,100];
s1 = SuperQuadrics({1+10*rand(3,1), 0.1+1.8*rand(2,1), -0.5+rand(1,2),...
    zeros(3,1), [1,0,0,0], mm});
s2 = SuperQuadrics({1+5*rand(3,1), 0.1+1.8*rand(2,1), -0.5+rand(1,2),...
    zeros(3,1), [1,0,0,0], mm});

% Transformations
M1 = quat2rotm(rand(1,4));
M2 = quat2rotm(rand(1,4));

% Compute Minkowski sums
mink_obj = MinkSumClosedForm(s1, s2, M1, M2);
m1 = s1.GetGradientsCanonical();
mink = mink_obj.GetMinkSumFromGradient(m1);

%% Plot: Minkowski sums
figure; hold on; axis equal; axis off;
X = reshape(mink(1,:), mm);
Y = reshape(mink(2,:), mm);
Z = reshape(mink(3,:), mm);
surf(X, Y, Z, 'FaceColor', 'y',...
    'EdgeColor', 'none', 'FaceAlpha', 0.6);

% Plot s1
x1 = M1 * s1.GetPointsCanonical();
X = reshape(x1(1,:), mm);
Y = reshape(x1(2,:), mm);
Z = reshape(x1(3,:), mm);
surf(X, Y, Z, 'FaceColor', 'g',...
    'EdgeColor', 'none', 'FaceAlpha', 1);

% Plot s2
for i = 1:size(mink,2)/5:size(mink,2)
    s2_touch = M2 * s2.GetPointsCanonical() + mink(:,i);
    
    X = reshape(s2_touch(1,:), mm);
    Y = reshape(s2_touch(2,:), mm);
    Z = reshape(s2_touch(3,:), mm);
    surf(X, Y, Z, 'FaceColor', 'b',...
        'EdgeColor', 'none', 'FaceAlpha', 1);
end

lightangle(gca,45,30);
lighting gouraud;