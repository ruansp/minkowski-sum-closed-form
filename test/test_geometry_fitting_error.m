% Test script for fitting approximation errors between different geometric
% models
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    convhull, sq_fitting, SuperQuadrics

close all; clear; clc;
add_paths();

N_trial = 100;

disp('****************************************************************')
disp('*** Testing for geometric model fitting approximation errors ***')
disp('****************************************************************')

%% Fitting superquadrics to convex polyhedra
disp('Fitting superquadrics to convex polyhedra...')
N_pts = 100;

vol_poly = nan(1, N_trial);
vol_sq = nan(1, N_trial);
rel_vol_sq_poly = nan(1, N_trial);
error_rel_sq_poly = nan(1, N_trial);
error_pt_surf = nan(1, N_trial);

for i = 1:N_trial
    % Generate random convex polyhedra
    point_cloud = [20;15;10] .* (2*rand(3,N_pts) - 1);
    [k, vol_poly(i)] = convhull(point_cloud(1,:), point_cloud(2,:),...
        point_cloud(3,:));
    poly = point_cloud(:,unique(k));
    
    % Fit superquadric model
    [sq_shape, error_pt_surf(i)] = sq_fitting([], poly);
    sq = SuperQuadrics({sq_shape.a, sq_shape.eps, [0,0],...
        sq_shape.pose.tc, sq_shape.pose.q, [50, 50]});
    vol_sq(i) = sq.GetVolume();
    
    poly_canonical = quat2rotm(sq.q)\(poly - sq.tc);
    error_pt_surf(i) = sum( abs(sq.GetImplicitFunction(poly_canonical)) )...
        / size(poly_canonical, 2);
    
    % Relative volume
    error_rel_sq_poly(i) = abs(vol_sq(i)-vol_poly(i))./vol_poly(i);
    rel_vol_sq_poly(i) = vol_sq(i)./vol_poly(i);
end

disp(['Mean error avg point to fitted surface distance: ',...
    num2str( mean(error_pt_surf) )]);
disp(['Mean error rel volume SQ fit Polyhedra: ',...
    num2str( mean(error_rel_sq_poly) )]);
disp(['Mean relative volume SQ fit Polyhedra: ',...
    num2str( mean(rel_vol_sq_poly) )]);

% Plots
figure; hold on; axis equal; axis off;
% plot3(point_cloud(1,:), point_cloud(2,:), point_cloud(3,:), 'bo');

plot3(poly(1,:), poly(2,:), poly(3,:), 'r*');
trisurf(k, point_cloud(1,:), point_cloud(2,:), point_cloud(3,:),...
    'FaceColor', 'y', 'FaceAlpha', 0.7)

sq.PlotShape('g', 0.3);

lightangle(gca,45,30);
lighting gouraud;

%% Fitting convex polyhedra to superquadrics
disp('****************************************************************')
disp('Fitting convex polyhedra to superquadrics...')
N_vtx = 4:20;

vol_sq_2 = nan(length(N_vtx), N_trial);
vol_poly_2 = nan(length(N_vtx), N_trial);
error_rel_poly_sq = nan(length(N_vtx), N_trial);

for m = 1:length(N_vtx)
    for i = 1:N_trial
        % Generate SQ
        sq2 = SuperQuadrics({[20,15,10].*rand(1,3), 2*rand(1,2), [0,0],...
            [0;0;0], rand(1,4), [N_vtx(m), N_vtx(m)]});
        vol_sq_2(m, i) = sq2.GetVolume();
        
        % Polyhedral approximation
        pts2 = sq2.GetPoints();
        [k, vol_poly_2(m, i)] = convhull(pts2');
        poly2 = pts2(:,unique(k));
        
        % Relative volume
        error_rel_poly_sq(m, i) = abs(vol_poly_2(m, i)-vol_sq_2(m, i))./...
            vol_sq_2(m,i);
    end
end

disp(['Mean rel volume Polyhedra fit SQ: ',...
    num2str( mean(error_rel_poly_sq, 2)' )]);

% Plots
% Fitting example
figure; hold on; axis equal; axis off;
sq2.N = [50, 50];
[sq2.eta, sq2.omega] = meshgrid(...
    -pi/2:pi/(sq2.N(1)-1):pi/2,...
    -pi-1e-6:2*pi/(sq2.N(2)-1):pi+1e-6);
sq2.PlotShape('g', 0.3);
plot3(pts2(1,:), pts2(2,:), pts2(3,:), 'bo');

plot3(poly2(1,:), poly2(2,:), poly2(3,:), 'r*');
trisurf(k, pts2(1,:), pts2(2,:), pts2(3,:),...
    'FaceColor', 'cyan', 'FaceAlpha', 0.3);

lightangle(gca,45,30);
lighting gouraud;

% Error metric of relative volume
figure; hold on; grid on;
boxplot(error_rel_poly_sq' * 100);

numVtx = N_vtx.^2;
for i = 1:size(numVtx,2)
    numVtxCell{i} = num2str(numVtx(i));
end

xtickangle(45)
xticklabels(numVtxCell)
% title('Relative volume for 100 random configurations')
xlabel('Number of vertices on surface')
ylabel('Relative volume error/%', 'Interpreter', 'tex')
set(gca,'FontSize', 16)

%% Fitting superquadrics and convex polyhedra to basic primitives
disp('****************************************************************')
disp('Fitting superquadrics and convex polyhedra to basic primitives...')

% Cube
disp('** Cube **')
vol_cube = nan(1, N_trial);
vol_poly_cube = nan(1, N_trial);
vol_sq_cube = nan(1, N_trial);
error_rel_sq_cube = nan(1, N_trial);
for i = 1:N_trial
    % Define cube
    cube_size = [50,36,20] .* rand(1,3) + 1;
    vol_cube(i) = prod(cube_size);
    
    % Polyhedron (exact)
    vol_poly_cube(i) = vol_cube(i);
    
    % SQ approximation
    sq_cube = SuperQuadrics({cube_size/2, [0.1,0.1], [0,0],...
            [0;0;0], [0,1,0,0], [50, 50]});
    vol_sq_cube(i) = sq_cube.GetVolume();
    
    % Relative volume
    error_rel_sq_cube(i) = abs(vol_sq_cube(i)-vol_cube(i))./...
        vol_cube(i);
end

disp(['Mean rel volume SQ fit Cube: ',...
    num2str( mean(error_rel_sq_cube) )]);

figure; hold on; axis equal;
sq_cube.PlotShape('g', 0.3);

lightangle(gca,45,30);
lighting gouraud;

% Ellipsoid
disp('** Ellipsoid **')
vol_ellip = nan(1, N_trial);
vol_poly_ellip = nan(1, N_trial);
vol_sq_ellip = nan(1, N_trial);
vol_rel_poly_ellip = nan(1, N_trial);
for i = 1:N_trial
    % Define ellipsoid
    ellips_size = [50,36,20] .* rand(1,3) + 1;
    vol_ellip(i) = 4/3 * pi * prod(ellips_size);
    
    % SQ (exact)
    sq_ellip = Ellipsoid({ellips_size, [0,0],...
            [0;0;0], [0,1,0,0], [9, 9]});
    vol_sq_ellip(i) = sq_ellip.GetVolume();
    
    % Polyhedral approximation
    pts_ellip = sq_ellip.GetPoints();
    [k, vol_poly_ellip(i)] = convhull(pts_ellip');
    
    % Relative volume
    vol_rel_poly_ellip(i) = abs(vol_poly_ellip(i)-vol_ellip(i))./...
        vol_ellip(i);
end

disp(['Mean rel volume Polyhedra fit Ellipsoid: ',...
    num2str( mean(vol_rel_poly_ellip) )]);

figure; hold on; axis equal; axis off;
plot3(pts_ellip(1,:), pts_ellip(2,:), pts_ellip(3,:), 'r*');
trisurf(k, pts_ellip(1,:), pts_ellip(2,:), pts_ellip(3,:),...
    'FaceColor', 'cyan', 'FaceAlpha', 0.3)

sq_ellip.N = [50, 50];
[sq_ellip.eta, sq_ellip.omega] = meshgrid(...
    -pi/2:pi/(sq_ellip.N(1)-1):pi/2,...
    -pi-1e-6:2*pi/(sq_ellip.N(2)-1):pi+1e-6);
sq_ellip.PlotShape('g', 0.3);

lightangle(gca,45,30);
lighting gouraud;

% Cylinder
disp('** Cylinder **')
N_poly_pts = 60;
vol_cyl = nan(1, N_trial);
vol_poly_cyl = nan(1, N_trial);
vol_sq_cyl = nan(1, N_trial);
vol_rel_poly_cyl = nan(1, N_trial);
vol_rel_sq_cyl = nan(1, N_trial);
for i = 1:N_trial
    % Define cylinder
    cyl_size = [50,36,20].*rand(1,3) + 1;
    vol_cyl(i) = pi * prod(cyl_size);
    
    % SQ approximation
    sq_cyl = SuperQuadrics({[cyl_size(1:2), cyl_size(3)/2], [0.1,1],...
        [0,0], [0;0;0], [0,1,0,0], [50, 50]});
    vol_sq_cyl(i) = sq_cyl.GetVolume();
    
    % Polyhedral approximation
    th = -pi:2*pi/(N_poly_pts/2-1):pi;
    cross_section = cyl_size(1:2)'.*[cos(th);sin(th)];
    pts_cyl = [[cross_section; -cyl_size(3)/2 *...
        ones(1,size(cross_section,2))], [cross_section;...
        cyl_size(3)/2 * ones(1,size(cross_section,2))]];
    [k, vol_poly_cyl(i)] = convhull(pts_cyl');
    
    % Relative volume
    vol_rel_poly_cyl(i) = abs(vol_poly_cyl(i)-vol_cyl(i))./...
        vol_cyl(i);
    vol_rel_sq_cyl(i) = abs(vol_sq_cyl(i)-vol_cyl(i))./...
        vol_cyl(i);
end

disp(['Mean rel volume Polyhedra fit Cylinder: ',...
    num2str( mean(vol_rel_poly_cyl) )]);
disp(['Mean rel volume SQ fit Cylinder: ',...
    num2str( mean(vol_rel_sq_cyl) )]);

figure; hold on; axis equal; axis off;
plot3(pts_cyl(1,:), pts_cyl(2,:), pts_cyl(3,:), 'r*');
trisurf(k, pts_cyl(1,:), pts_cyl(2,:), pts_cyl(3,:),...
    'FaceColor', 'cyan', 'FaceAlpha', 0.3)

sq_cyl.PlotShape('g', 0.3);

lightangle(gca,45,30);
lighting gouraud;
