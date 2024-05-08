%% Main 04 -- Manual implementation of Full Horizon Estimation.
%
% The goal of this script is to show:
%
% 1. How Full Horizon Estimation works.
% 2. How to implement manually Full Horizon Estimation.

close all;
clear all;

%% 1. How Full Horizon Estimation works.
%
% In the previous script, 'main_03', we found the initial condition that best
% explains the experimental dataset. 
%
% However, to do that, we used all the measurements of the experiment.
% Therefore, using that method, we have to wait until the experiement has ended
% (until we have all the measurements) to get the estimates of the states. This
% is an off-line estimation method.
%
% If we want to have on-line estimantion, we have to use a different approach
% to define the windows we want to solve. This is the orgin of the full-horizon
% estimation: an on-line estimation method which uses all the available
% measurements.
%
% The idea is a follows:
%
% 1. At time 'k', define a window that covers all measurements from 'k0' (the
%    beggining of the experiment) up to current time ('k').
% 2. Solve the window with the available measurements.
% 3. At time 'k+1', define a new window that covers also that new available 
%    measurement.
% 4. Solve the window 'k+1'.
% 5. Repeat until this procedure until the expriment finish.

%% 2. How to implement manually Full Horizon Estimation.
%
% Now we will implement full-horizon estimation using the previous code.

% First we load the estimation problem as we did in 'main03'. However, this
% time we will use a helper function, 'load_estimation_problem', to do this.
data = "exp01";
problem = utils.load_estimation_problem(data);

% Second, we will plot the dataset using also a helper function.
f = figure();
utils.plot_data(data);

% Define and solve iteratively all the windows for Full Horizon Estimation.
for k = 1:length(problem.t)
    % In FHE, 'k0' is always the first measurement.
    k0 = 1;

    % A priori initial condition of the window 'k'.
    x0_k = problem.x0;
    
    % Covariance of the a priori initial condition of the window 'k'.
    Q0_k = problem.Q0;
 
    % Definition of window 'k'.
    window = mohotools.Window(problem, [k0 k], x0_k, Q0_k);

    % Solve window 'k'.
    solution = window.solve();
    
    % Plot window 'k'.
    figure(f);
    utils.plot_window(solution, ['FHE-' num2str(k)]);
end
