% Demo script for exact closed-form Minkowski sums between two 3D 
% superquadrics that represents basic geometric primitives
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    defineSQFromPrimitives, MinkSumClosedForm

close all; clear; clc;
add_paths();
mm = [50,50];

disp('******************************************************************')
disp('* Demonstration of closed-form Minkowski sums between primitives *')
disp('******************************************************************')

%% Geometric primitives
% Box
box1 = [5.2, 3.3, 2.87];
box2 = [4.3, 2.8, 1.254];

% Cylinder
cyl1 = [6.24, 3.6, 10];
cyl2 = [3.58, 2.52, 5];

% Sphere
sph1 = 8.6 * ones(1,3);
sph2 = 3.56 * ones(1,3);

% Parallelepiped
M1 = [1, 0, 0.5; 0, 1, 0.2; 0, 0, 1];
M2 = [1, 0, 0.47; 0, 1, 0.14; 0, 0, 1];
% M2 = eye(3);

%% Superquadric models
% -> Please uncomment the lines for primitives that you would like to see

s1 = defineSQFromPrimitive('box', box1, mm);
% s1 = defineSQFromPrimitive('cylinder', cyl1, mm);
% s1 = defineSQFromPrimitive('ellipsoid', sph1, mm);

s2 = defineSQFromPrimitive('box', box2, mm);
% s2 = defineSQFromPrimitive('cylinder', cyl2, mm);
% s2 = defineSQFromPrimitive('ellipsoid', sph2, mm);

%% Compute Minkowski sums
mink_obj_shear = MinkSumClosedForm(s1, s2,...
    quat2rotm(s1.q), quat2rotm(s2.q));
m1 = s1.GetGradientsCanonical();
mink = mink_obj_shear.GetMinkSumFromGradient(m1);

% Plots
figure; hold on; axis equal; axis off;
% Minkowski sums
X = reshape(mink(1,:), mm);
Y = reshape(mink(2,:), mm);
Z = reshape(mink(3,:), mm);
surf(X, Y, Z, 'FaceColor', 'y',...
    'EdgeColor', 'none', 'FaceAlpha', 0.6);

% s1 and s2
s1.PlotShape('g', 1);
s2.tc = mink(:,1000);
s2.PlotShape('b', 1);

lightangle(gca,45,30);
lighting gouraud;

%% Linear transform of original bodies
mink_obj_shear = MinkSumClosedForm(s1, s2,...
    quat2rotm(s1.q) * M1, quat2rotm(s2.q) * M2);
m1 = s1.GetGradientsCanonical();
mink = mink_obj_shear.GetMinkSumFromGradient(m1);

% Plots
figure; hold on; axis equal; axis off;
% Minkowski sums
X = reshape(mink(1,:), mm);
Y = reshape(mink(2,:), mm);
Z = reshape(mink(3,:), mm);
surf(X, Y, Z, 'FaceColor', 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.6);

% s1 and s2
x1_shear = quat2rotm(s1.q) * M1 * s1.GetPointsCanonical();
x = reshape(x1_shear(1,:), mm);  
y = reshape(x1_shear(2,:), mm);
z = reshape(x1_shear(3,:), mm);
surf(x, y, z, 'FaceColor', 'g', 'FaceAlpha', 1, 'EdgeColor', 'none');

x2_shear = quat2rotm(s2.q) * M2 * s2.GetPointsCanonical() + mink(:,1000);
x = reshape(x2_shear(1,:), mm);  
y = reshape(x2_shear(2,:), mm);
z = reshape(x2_shear(3,:), mm);
surf(x, y, z, 'FaceColor', 'b', 'FaceAlpha', 1, 'EdgeColor', 'none');

lightangle(gca,45,30);
lighting gouraud;

%% Define superquadric model
function sq = defineSQFromPrimitive(geom_type, geom_size, mm)
% defineSQFromPrimitive Define a 3D superquadric model for geometric
% primitives, e.g. box, cylinder and ellpsoid

sq = SuperQuadrics({[1,1,1], [1,1], [0,0],...
    zeros(3,1), rand(1,4), mm});

switch geom_type
    case 'box'
        sq.a = geom_size/2;
        sq.eps = [0.1, 0.1];
    case 'cylinder'
        sq.a = [geom_size(1:2), geom_size(3)/2];
        sq.eps = [0.1, 1];
    case 'ellipsoid'
        sq.a = geom_size;
        sq.eps = [1, 1];
end

disp(['** Geometry type: ', geom_type, ' ** size: ', num2str(geom_size)])

end