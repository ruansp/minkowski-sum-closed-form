% Benchmark file for closed-form Minkowski sums between two 2D SQs
%  Algorithms compared: convex hull, edge sorting, closed-form geometric,
%  closed-form normal/gradient (proposed)
%
%  Author: 
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    SuperEllipse, Ellipse, MinkSumClosedForm, MinkSumDefinition

close all; clear; clc;
add_paths();

num_vtx = 20:10:500;
num_trials = 100;

time_closed = nan(num_trials, length(num_vtx));
time_closed_geom = nan(num_trials, length(num_vtx));
time_edge_sort = nan(num_trials, length(num_vtx));
time_convhull = nan(num_trials, length(num_vtx));

%% Geometric type
% -> Please uncomment for different geometries

geom_type = 'Ellipse-Ellipse';
% geom_type = 'Ellipse-Superellipse';
% geom_type = 'Superellipse-Superellipse';

disp(['Minkowski sum geometric type: ', geom_type]);

%% Benchmark
for mm = 1:length(num_vtx)
    disp(['Num of vertices: ', num2str(num_vtx(mm))]);
    
    for nn = 1:num_trials
        s1 = SuperEllipse([1+10*rand(1,2),0.1+1.8*rand,0,0,0,0,num_vtx(mm)]);
        s2 = SuperEllipse([1+10*rand(1,2),0.1+1.8*rand,0,0,0,0,num_vtx(mm)]);
        
        % Transformations for the two SQs
        angles = 2*pi*rand(1,2);
        shears = 2*(2*rand(1,2)-1);
        
        M1 = rot2(angles(1)) * [1, shears(1); 0, 1];
        M2 = rot2(angles(2)) * [1, shears(2); 0, 1];
        
        if strcmp(geom_type, 'Ellipse-Ellipse')
            s1 = Ellipse([1+10*rand(1,2),0,0,0,0,num_vtx(mm)]);
            s2 = Ellipse([1+10*rand(1,2),0,0,0,0,num_vtx(mm)]);
            
        elseif strcmp(geom_type, 'Ellipse-Superellipse')
            s2 = Ellipse([1+10*rand(1,2),0,0,0,0,num_vtx(mm)]);
        end
        
        %% Closed-form Minkowski sums (proposed)
        minkObj = MinkSumClosedForm(s1, s2, M1, M2);
        
        t_closed= tic;
        m1 = s1.GetGradients();
        mink = minkObj.GetMinkSumFromGradient(m1);
        time_closed(nn,mm) = toc(t_closed);
        
        %% Closed-form Minkowski sums using geometric method
        % Only works for at least one ellipsoid
        if ~strcmp(geom_type, 'Superellipse-Superellipse')
            t_closed_geom = tic;
            mink_geom = minkObj.GetMinkSumGeometric();
            time_closed_geom(nn,mm) = toc(t_closed_geom);
        end
        
        %% Minkowski sums using edge sorting algorithm
        t_edge_sort = tic;
        mink_edge_sort = MinkSumEdgeSort2D(s1, s2, M1, M2);
        time_edge_sort(nn,mm) = toc(t_edge_sort);
        
        %% Minkowski sums using definition      
        t_def = tic;
        mink_def = MinkSumDefinition(s1, s2, M1, M2);
        time_convhull(nn,mm) = toc(t_def);
    end
end

%% Plots
fontSize = 15;

% Time over number of vertices
figure; hold on; grid on;
plot(num_vtx, mean(time_closed,1), '-b', 'LineWidth', 2);
plot(num_vtx, mean(time_edge_sort,1), '-*g', 'LineWidth', 2);
plot(num_vtx, mean(time_convhull,1), '--c', 'LineWidth', 2);
plot(num_vtx, mean(time_closed_geom,1), '-.r', 'LineWidth', 2);

shadedErrorBar(num_vtx, mean(time_closed,1), std(time_closed,1),...
    'lineprops','-b')
shadedErrorBar(num_vtx, mean(time_edge_sort,1),...
    std(time_edge_sort,1),...
    'lineprops','-*g')
shadedErrorBar(num_vtx, mean(time_convhull,1), std(time_convhull,1),...
    'lineprops','--c')
shadedErrorBar(num_vtx, mean(time_closed_geom,1),...
    std(time_closed_geom,1),...
    'lineprops','-.r')

xlabel('Number of vertices')
ylabel('Running time/s')

if ~strcmp(geom_type, 'Superellipse-Superellipse')
    legend('Closed-form (proposed)', 'Edge sort', 'Convex hull',...
        'Closed-form Geometric')
else
    legend('Closed-form (proposed)', 'Edge sort', 'Convex hull')
end

set(gca,'fontsize',fontSize);