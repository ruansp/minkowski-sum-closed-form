% Verification script for the proposed closed-form Minkowski sums between 
% two 3D superquadrics
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    kiss_point_measure

close all; clear; clc;
add_paths();

% Parameters
num_vtx = [50,50];
num_trials = 100;

disp('Start verification...')
dist = nan(num_trials, 1);
num_contacts = nan(num_trials, 1);
err_kiss_point = nan(num_trials, 2);

for nn = 1:num_trials
    s1 = SuperQuadrics({1+9*rand(3,1), 0.1+1.8*rand(2,1), [0,0],...
        zeros(3,1), [1,0,0,0], num_vtx});
    s2 = SuperQuadrics({1+5*rand(3,1), 0.1+1.8*rand(2,1), [0,0],...
        zeros(3,1), [1,0,0,0], num_vtx});
    
    % Transformations
    M1 = quat2rotm(pi*rand(1,4));
    M2 = quat2rotm(pi*rand(1,4));
    
    %% Compute Minkowski sums
    minkObj = MinkSumClosedForm(s1, s2, M1, M2);
    x1 = M1 * s1.GetPoints();
    m1 = s1.GetGradients();
    mink = minkObj.GetMinkSumFromGradient(m1);
    
    %% Verification metrics
    % Closed-form boundary points and definition distance
    [num_contacts(nn), dist(nn)] = avg_dist_implicit(...
        s1, s2, M1, M2, mink);
    
    % Kissing point verification
    err_kiss_point(nn, :) = kiss_point_measure(s1, s2, M1, M2, x1, mink);
    
    clc
    disp(['Progress: ', num2str(nn/num_trials*100),'%'])
end

%% Results
disp('**************************************************************')
disp('*** Verifications of closed-form Minkowki sums for 3D case ***')
disp('**************************************************************')
disp('Averaged distance from implicit surface (DI):')
disp(['Mean: ', num2str(mean(dist,1))])
disp(['Standard deviation: ', num2str(std(dist,1))])

disp('************************************************')
disp('Contact number (Nc):')
disp(['Mean: ', num2str(mean(num_contacts,1))])
disp(['Standard deviation: ', num2str(std(num_contacts,1))])

disp('************************************************')
disp('Kissing point measure [implicit function, gradient]:')
disp(['Mean: ', num2str(mean(err_kiss_point,1))])
disp(['Standard deviation: ', num2str(std(err_kiss_point,1))])

%% Plots
fontSize = 15;

% Results based on distance from implicit surface
figure; hold on; grid on;
plot(1:num_trials, dist, '-*', 'LineWidth', 2)

xlabel('Index of trials')
ylabel('$D_I$', 'Interpreter', 'latex')

set(gca,'fontsize',fontSize);

% Results based on kissing point measure
figure; hold on; grid on;
plot(1:num_trials, err_kiss_point, '*', 'LineWidth', 2)

xlabel('Index of trials')
ylabel('$e_{kp}$', 'Interpreter', 'latex')
legend('Implicit', 'Gradient')

set(gca,'fontsize',fontSize);