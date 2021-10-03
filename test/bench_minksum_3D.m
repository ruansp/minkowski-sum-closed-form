% Benchmark script for Minkowski sums between two 3D SQs
%  Algorithms compared: convex hull, closed-form geometric, closed-form
%  normal/gradient (proposed)
%
%  Author: 
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    SuperQuadrics, Ellipsoid, MinkSumClosedForm, MinkSumDefinition

close all; clear; clc;
add_paths();

num_vtx = (5:5:50)' * ones(1,2);
num_trials = 100;

time_closed = nan(num_trials, length(num_vtx));
time_closed_geom = nan(num_trials, length(num_vtx));
time_convhull = nan(num_trials, length(num_vtx));

%% Geometric type
% -> Please uncomment for different geometries

geom_type = 'Ellipsoid-Ellipsoid';
% geom_type = 'Ellipsoid-Superquadric';
% geom_type = 'Superquadric-Superquadric';

disp(['Minkowski sum geometric type: ', geom_type]);

%% Benchmark
for mm = 1:length(num_vtx)
    disp(['Num of vertices: ', num2str( prod(num_vtx(mm,:)) )]);
    
    for nn = 1:num_trials
        % Define two body shapes
        s1 = SuperQuadrics({1+9*rand(3,1), 0.2+1.8*rand(2,1), [0,0],...
            zeros(3,1), [1,0,0,0], num_vtx(mm,:)});
        s2 = SuperQuadrics({1+5*rand(3,1), 0.2+1.8*rand(2,1), [0,0],...
            zeros(3,1), [1,0,0,0], num_vtx(mm,:)});
        
        if strcmp(geom_type, 'Ellipsoid-Ellipsoid')
            s1 = Ellipsoid({1+9*rand(3,1), [0,0],...
                zeros(3,1), [1,0,0,0], num_vtx(mm,:)});
            s2 = Ellipsoid({1+5*rand(3,1), [0,0],...
                zeros(3,1), [1,0,0,0], num_vtx(mm,:)});
            
        elseif strcmp(geom_type, 'Ellipsoid-Superquadric')
            s2 = Ellipsoid({1+5*rand(3,1), [0,0],...
                zeros(3,1), [1,0,0,0], num_vtx(mm,:)});
        end
        
        % Transformations
        M1 = quat2rotm(pi*rand(1,4));
        M2 = quat2rotm(pi*rand(1,4));
        
        %% Compute closed-form Minkowski sums
        minkObj = MinkSumClosedForm(s1, s2, M1, M2);
        
        t_closed = tic;
        m1 = s1.GetGradients();
        mink = minkObj.GetMinkSumFromGradient(m1);
        time_closed(nn,mm) = toc(t_closed);
        
        %% Compute closed-form Minkowski sums using geometric method
        % Only works for at least one ellipsoid
        if ~strcmp(geom_type, 'Superquadric-Superquadric')
            t_closed_geom = tic;
            mink_geom = minkObj.GetMinkSumGeometric();
            time_closed_geom(nn,mm) = toc(t_closed_geom);
        end
        
        %% Compute Minkowski sums using definition
        t_def = tic;
        mink_def = MinkSumDefinition(s1, s2, M1, M2);
        time_convhull(nn,mm) = toc(t_def);
    end
end

%% Plots
fontSize = 15;

% Time over number of vertices
figure; hold on; grid on;
plot(num_vtx(:,1).*num_vtx(:,2), mean(time_closed,1),...
    '-b', 'LineWidth', 2);
plot(num_vtx(:,1).*num_vtx(:,2), mean(time_convhull,1),...
    '--c', 'LineWidth', 2);
plot(num_vtx(:,1).*num_vtx(:,2), mean(time_closed_geom,1),...
    '-.r', 'LineWidth', 2);

shadedErrorBar(num_vtx(:,1).*num_vtx(:,2), mean(time_closed,1),...
    std(time_closed,1),...
    'lineprops','-b')
shadedErrorBar(num_vtx(:,1).*num_vtx(:,2), mean(time_convhull,1),...
    std(time_convhull,1),...
    'lineprops','--c')
shadedErrorBar(num_vtx(:,1).*num_vtx(:,2), mean(time_closed_geom,1),...
    std(time_closed_geom,1),...
    'lineprops','-.r')

xlabel('Number of vertices')
ylabel('Running time/s')

if ~strcmp(geom_type, 'Superquadric-Superquadric')
    legend('Closed-form (proposed)', 'Convex hull', 'Closed-form Geometric')
else
    legend('Closed-form (proposed)', 'Convex hull')
end

set(gca,'fontsize',fontSize);