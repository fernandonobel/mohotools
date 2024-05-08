%% Main 05 -- Manual implementation of Moving Horizon Estimation.
%
% The goal of this script is to show:
%
% 1. How Moving Horizon Estimation works.
% 2. How to implement manually Moving Horizon Estimation.

close all;
clear all;

%% 1. How Moving Horizon Estimation works.
%
% Full-horizon estimations allows us to have on-line estimations. However, the
% size of the window we are solving increases with the number of available
% measurements. The problem is that the time it takes to solve a window is
% directly realted to its size. Therefore, in experiments where there are
% available a lot of experimental measurements, full horizon estimation gets
% progressibely slower up to the point where it not longer possible to apply
% this method on-line.
%
% The solution to this problem is moving-horizon estimation. In this case, we
% limit the maximum size of a window to 'N' (a window of size 'N' covers only
% the 'N+1' most recent measurements). This way we avoid the unbound growth of
% the window and we ensure that the method can be applied on-line.
%
%
% The idea is:
%
% - When the number of available measurements is less than or equal to 'N+1', 
%   run full-horizon estimation.
% - When the number of available measurements is larger than 'N+1', limit the
%   window size to 'N' (run moving-horizon estimation).

%% 2. How to implement manually Moving Horizon Estimation.
%
% The implementation of moving-horizon estimation is equal to the full-horizon
% version, but we need to limit the window size.

%% Load the estimation problem.

data = "exp01";

f = figure();
utils.plot_data(data);

problem = utils.load_estimation_problem(data);

% Window size.
N = 5;

%% Run Full Horizon Estimation (FHE) at the beggining.

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
    solution = window.solve();
    
    % Plot window 'k'.
    figure(f);
    utils.plot_window(solution, ['FHE-' num2str(k)]);
end

%% Run Moving Horizon Estimation (MHE) for the rest of the experiment.

for k = N+2 : length(problem.t)
    
    % In MHE, 'k0' is 'N' measurements behind 'k'. 
    k0 = k-N;

    % We don't have the a prori initial condition for time 'k'.
    % Let's use the a priori initial condition of the problem definition.
    x_k0 = problem.x0;
    
    % But, increase its covariance as it is wrong.
    Q_k0 = problem.Q0*5;
    
    % Time vector for simulating the solution.
    tspan = [k0:k+1] - 1;
    
    % Definition of window 'k'.
    window = mohotools.Window(problem, [k0 k], x_k0, Q_k0, tspan);
    
    % Solve window 'k'.
    solution = window.solve();
    
    % Plot window 'k'.
    figure(f);
    utils.plot_window(solution, ['MHE-' num2str(k)]);
end
