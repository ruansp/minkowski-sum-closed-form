% Main script for application of c-obstacle generations using closed-form
% Minkowski sums, SE(3)
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    SuperQuadrics, generate_obstacle, MinkSumClosedForm
%
%  External link
%    SO(3) sampling using double-coset decomposition
%    (https://github.com/ruansp/quantization-double-coset)

close all; clear; clc;

add_paths();

% Load SO(3) samples from double-coset decomposition
load('../misc/R_coset.mat');

% Robot
robot = SuperQuadrics({5*rand(3,1), 2*rand(2,1), [0,0], zeros(3,1),...
    [1,0,0,0], [10,10]});
robot_config = [robot.a', robot.eps', robot.tc',...
    robot.q];

% Generate environment
N = [50,200,1e3];
N_rot = size(R_coset,3);
c_obs = cell(1,length(N));
obstacle = cell(1,length(N));
obs_config = cell(1,length(N));

run_time = zeros(1,length(N));

%% Main routine
disp('****************************************************************')
disp('Application on configuration-space obstacle generations in SE(3)')
disp('****************************************************************')
for i = 1:length(N)
    [obstacle{i}, obs_config{i}] = generate_obstacle(N(i));
    
    % Minkowski sums
    mink = cell(N(i), N_rot);
    ti = tic;
    for k = 1:N_rot
        if mod(k/N_rot*100,1) == 0
            clc
            disp(['Number of obstacles: ', num2str(N(i)),...
                ' -- Progress: ', num2str(k/N_rot*100), '%'])
        end
        
        for j = 1:N(i)
            minkObj = MinkSumClosedForm(robot,obstacle{i}{j},...
                R_coset(:,:,k), quat2rotm(obstacle{i}{j}.q));
            m1 = robot.GetGradients();
            mink{j,k} = minkObj.GetMinkSumFromGradient(m1) +...
                obstacle{i}{j}.tc;
        end
    end
    run_time(i) = toc(ti);
    
    c_obs{i} = mink;
end

%% Plots
figure; hold on; axis off;
idx = 3;
k_rot = 15;
for i = 1:N(idx)
    obstacle = SuperQuadrics({obs_config{idx}(i,1:3)',...
        obs_config{idx}(i,4:5)', [0,0], obs_config{idx}(i,6:8)',...
        obs_config{idx}(i,9:12), [10,10]});
    obstacle.PlotShape('b', 0.7);
    
    plot3(c_obs{idx}{i,k_rot}(1,:), c_obs{idx}{i,k_rot}(2,:),...
        c_obs{idx}{i,k_rot}(3,:), 'r.', 'LineWidth', 0.1)
end

figure; hold on;
plot(N, run_time);

%% Subroutine for generate N random obstacles
function [obstacle, obs_config] = generate_obstacle(N)
obs_a = 10 * rand(3,N);
obs_eps = 2 * rand(2,N);
obs_tc = 200 * (2*rand(3,N)-1);
obs_quat = pi * rand(4,N);

obstacle = cell(1,N);
for i = 1:N
    obstacle{i} = SuperQuadrics({obs_a(:,i), obs_eps(:,i), [0,0],...
        obs_tc(:,i), obs_quat(:,i)', [10,10]});
end

obs_config = [obs_a; obs_eps; obs_tc; obs_quat]';
end