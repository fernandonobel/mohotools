%% Main 06 -- Arrival cost
%
% The goal of this script is to show:
%
% 1. How the arrival cost improves the performance of MHE.
% 2. How to implment maximum-likelihood arrival cost.

close all;
clear all;

%% 1. How the arrival cost improves the performance of MHE.
%
% The window, in moving-horizon estimation, is limited in size. Therefore, the 
% window does not use all the avilable measurements (it only used the N+1 most 
% recent measurements). This means than we are loosing past information and
% this causes a loss of performance in moving-horizon estimation with respect
% full-horizon estimation.
%
% One way to fix this problem is calculating the arrival cost. The arrival cost
% is the best guess we can get for the a priori initial condition of a window
% using only the measurements that are outside that window. This way, the
% arrival cost is a way to re-introduce past information in the current window.
%
% How can we calculate the arrival cost? 
%
% Typically, the arrival cost is not calculated (as it would require to solve
% the full-horizon estimation problem). Instead, we _approximate_ the arrival
% cost.
%
% There are many different methods to approximate the arrival cost. In this
% script, we show the maximum-likelihood method. The arrival cost for window
% 'k' is approximated from the solution of the previous window 'k-N-1'. Becasue
% the estimation for the state at time 'k-N' of the window 'k-N-1' is a summary
% of all the past information outside of the current window 'k'. This method is
% called "maximum-likelihood" because the covariance of the arrival cost is
% directly the covariance of the state estimation using maximum-likelihood.

%% 2. How to implment maximum-likelihood arrival cost.
%
% To implement the maximum-likelihood arrival cost we have to store the
% solutions of previous windows and extend the simulation of the solution of 
% each window up to 'k+1'.

%% Load the estimation problem.

data = "exp01";

f = figure();
utils.plot_data(data);

problem = utils.load_estimation_problem(data);

% Window size.
N = 5;

% Cell array of the solutions of each window.
solution = {};

%% Full Horizon Estimation (FHE).

for k = 1 : N+1
    
    % In FHE, 'k0' is always the first measurement.
    k0 = 1;
    
    % A priori initial condition of the window 'k'.
    x_k0 = problem.x0;
    
    % Covariance of the a priori initial condition of the window 'k'.
    Q_k0 = problem.Q0;
    
    % Time vector for simulating the solution.
    tspan = [k0:k+1] - 1;
    
    % Definition of window 'k'.
    window = mohotools.Window(problem, [k0 k], x_k0, Q_k0, tspan);
    
    % Solve window 'k'.
    solution{k} = window.solve();
    
    % Plot window 'k'.
    figure(f);
    utils.plot_window(solution{k}, ['FHE-' num2str(k)]);
end

%% Moving Horizon Estimation (MHE).

for k = N+2 : length(problem.t)
    
    % In MHE, 'k0' is N measurements behind 'k'. 
    k0 = k-N;
    
    % A priori initial condition of the window 'k'.
    x_k0 = solution{k-N-1}.x(end,:)';
    
    % Covariance of the a priori initial condition of the window 'k'.
    Q_k0 = squeeze(solution{k-N-1}.P(end,:,:));
    
    % Time vector for simulating the solution.
    tspan = [k0:k+1] - 1;
    
    % Definition of window 'k'.
    window = mohotools.Window(problem, [k0 k], x_k0, Q_k0, tspan);
    
    % Solve window 'k'.
    solution{k} = window.solve();
    
    % Plot window 'k'.
    figure(f);
    utils.plot_window(solution{k}, ['MHE-' num2str(k)]);
end
