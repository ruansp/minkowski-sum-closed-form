% Main script for application on collision detection between two 
% superquadrics and benchmarks with other algorithms
%
%  Algorithms compared:
%    1. MinkSum-based (proposed)
%    2. Implicit surface, interior-point method
%    3. Parametric surface, common-normal concept
%    4. GJK (approximated faceted surface)
%    5. Algebraic separation condition (for ellipsoids only)
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    SuperQuadrics, Ellipsoid, collision_minksum, 
%    collision_implicit_interior, collision_parametric, GJK, 
%    collision_ellipsoid_asc

close all; clear; clc;

add_paths();
num_trials = 1e4;

time_mink_ray = nan(1,num_trials);
time_mink_normal = nan(1,num_trials);
time_implicit = nan(1,num_trials);
time_param = nan(1,num_trials);
time_gjk = nan(1,num_trials);
time_asc = nan(1,num_trials);

status_mink_ray = nan(1,num_trials);
status_mink_normal = nan(1,num_trials);
status_implicit = nan(1,num_trials);
status_param = nan(1,num_trials);
status_gjk = nan(1,num_trials);
status_asc = nan(1,num_trials);

pt_mink_ray = cell(1,num_trials);
pt_mink_normal = cell(1,num_trials);
pt_implicit = cell(1,num_trials);
pt_param = cell(1,num_trials);

%% Geometric types
% -> Please uncomment for different geometric types

% geom_type = "Ellipsoid";
geom_type = "SuperQuadrics";

if strcmp(geom_type, 'SuperQuadrics')
    algorithm_names = {'Mink (Ray)', 'Mink (Normal)', 'Implicit', ...
        'Parametric', 'GJK'};
elseif strcmp(geom_type, 'Ellipsoid')
    algorithm_names = {'Mink (Ray)', 'Mink (Normal)', 'Implicit', ...
        'Parametric', 'GJK', 'ASC'};
end

%% Benchmark
disp('****************************************************************')
disp('*** Application on collision detection between superquadrics ***')
disp('****************************************************************')

mm = [10,10];
figure; hold on; axis equal;
for i = 1:num_trials
    if mod(i/num_trials*100,1) == 0
        clc
        disp(['Progress: ', num2str(i/num_trials*100),'%'])
    end
    
    if strcmp(geom_type, 'SuperQuadrics')
        s1 = SuperQuadrics({1+10*rand(3,1), 0.1+1.8*rand(2,1), [0,0],...
            5*(2*rand(3,1)-1), rand(1,4), mm});
        s2 = SuperQuadrics({1+5*rand(3,1), 0.1+1.8*rand(2,1), [0,0],...
            10+5*(2*rand(3,1)-1), rand(1,4), mm});
    elseif strcmp(geom_type, 'Ellipsoid')
        s1 = Ellipsoid({1+10*rand(3,1), [0,0],...
            5*(2*rand(3,1)-1), rand(1,4), mm});
        s2 = Ellipsoid({1+5*rand(3,1), [0,0],...
            10+5*(2*rand(3,1)-1), rand(1,4), mm});
    end
    
    s1_surf = s1.GetSurf();
    s2_surf = s2.GetSurf();
    
    %% Minkowski-based (ray) status
    t_s = tic;
    [status_mink_ray(i), dist_mink_ray(i), pt_mink_ray{i},...
        f_mink_ray(i,:)] = collision_minksum(s1, s2, 'ray');
    time_mink_ray(i) = toc(t_s);
    
    %% Minkowski-based (common normal) status
    t_s = tic;
    [status_mink_normal(i), dist_mink_normal(i), pt_mink_normal{i},...
        f_mink_normal(i,:)] = collision_minksum(s1, s2, 'common-normal');
    time_mink_normal(i) = toc(t_s);
    
    %% Implicit surface status
    t_s = tic;
    [status_implicit(i), dist_implicit(i), pt_implicit{i}] =...
        collision_implicit_interior(s1, s2);
    time_implicit(i) = toc(t_s);
    
    %% Param surface common normal status
    t_s = tic;
    [status_param(i), dist_param(i), pt_param{i}, f_param(i,:)] =...
        collision_parametric(s1, s2);
    time_param(i) = toc(t_s);
    
    %% GJK status
    t_s = tic;
    status_gjk(i) = GJK(s1_surf, s2_surf, 100);
    time_gjk(i) = toc(t_s);
    
    %% ASC for ellipsoid
    if strcmp(geom_type, 'Ellipsoid')
        t_s = tic;
        status_asc(i) = collision_ellipsoid_asc(s1, s2);
        time_asc(i) = toc(t_s);
    end
end

%% Plots
% Time comparisons
figure;
if strcmp(geom_type, 'SuperQuadrics')
    boxplot([time_mink_ray', time_mink_normal',...
        time_implicit', time_param', time_gjk'],...
        'Labels', algorithm_names, 'OutlierSize', 0.01);
elseif strcmp(geom_type, 'Ellipsoid')
    boxplot([time_mink_ray', time_mink_normal',...
        time_implicit', time_param', time_gjk', time_asc'],...
        'Labels', algorithm_names, 'OutlierSize', 0.01);
end
ylim([0,0.05])

% Status comparisons
figure; hold on;
plot(1:num_trials, status_mink_ray, 'b*')
plot(1:num_trials, status_mink_normal, 'b--')
plot(1:num_trials, status_implicit, 'rd')
plot(1:num_trials, status_param, 'g.-')
plot(1:num_trials, status_gjk, 'ko')
plot(1:num_trials, status_asc, 'yo')
legend(algorithm_names)

%% Display: Objective success rate
func_thres = 1e-10;

func_mink_ray = sqrt(sum(f_mink_ray.^2, 2));
ratio_func_mink_ray = sum(func_mink_ray < func_thres) / num_trials

func_mink_normal = sqrt(sum(f_mink_normal.^2, 2));
ratio_func_mink_normal = sum(func_mink_normal < func_thres) / num_trials

func_param = sqrt(sum(f_param.^2, 2));
ratio_func_param = sum(func_param < func_thres) / num_trials

%% Display: Mean running time
mean_time_mink_ray = mean(time_mink_ray)
mean_time_mink_normal = mean(time_mink_normal)
mean_time_implicit = mean(time_implicit)
mean_time_param = mean(time_param)
mean_time_gjk = mean(time_gjk)
mean_time_asc = mean(time_asc)