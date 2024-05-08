%% Main 03 -- Solve one window.
%
% The goal of this script is to show:
%
% 1. How an estimation problem is defined in Mohotools.
% 2. How can we solve a window with Mohotools.

close all;
clear all;

%% 1. How the estimation problem is defined in Mohotools.
%
% Mohotools separates the definition of the problem we want to solve (which is
% the estimation of the value of a unknown state) from the method we want to
% use to solve it (for example, the full- or moving-horizon estimation methods).
%
% Mohotools defines a class, EstimationProblem, for handling all the
% information needed to define and solve a estimation problem. In this section,
% we will so how to create the EstimationProblem and how to fill it.
%
% The estimation problem is composed of four elements:
%
% a. The experimental measurements.
% b. The a priori knowledge of the experiment's initial condition.
% c. The state vector bounds.
% d. The model ode and output functions.

% We have to create an object of EstimationProblem class.
problem = mohotools.EstimationProblem();

%% 1.a Load the experimental measurements and the covariance matrix.

% Load the data.
measurements = readtable('data/exp01/measurements.csv');
measurements_covariance = readmatrix('data/exp01/measurements_covariance.csv');

% Experiment time vector.
problem.t = measurements.time;

% Experiment measurements.
problem.y = measurements.S;

% Covariance matrix of the measurements.
problem.R = measurements_covariance;

%% 1.b The a priori knowledge of the experiment's initial condition.

% A priori initial condition.
problem.x0 = [12; 0.5];

% Covariance matrix of the a priori initial condition.
problem.Q0 = diag([1 1]);

%% 1.c The state vector bounds.

% States upper_bound.
problem.x_ub = [20; 20];

% States lower bound.
problem.x_lb = [0; 0];

%% 1.d The model ode and output functions.

% Load parameter vector.
p_table = readtable("+model/parameter.csv");
p = p_table.value;

% Function handle which evaluates the model ode.
problem.ode = @(t,x) model.ode(t, x, p);

% Function handle which evaluates the model output.
problem.output = @(t,x) model.output(t, x, p);

%% 2. How can we solve a given window with Mohotools.
%
% The fundamental concept of full- and moving-horizon estimation is the concept
% of the window. The window is the set of measurements that we are considering
% when we estimate the intial condition of the experiment. 
%
% How we define the range of these windows (and how we connect their results) 
% leads to the formulation of the full- or moving-horizon estimation. However,
% as an introduction to implement full- and moving-horizon estimaton, let's
% define and solve one window.
%
% Next, we will define a window that covers the whole dataset and we will
% solve this window to find the optimal intitial condition.

%% 2.a Define the window.

% The window will start on the first measurement of the dataset.
k0 = 1;

% And the window will end on the last measurement.
k = length(problem.t);

% The a priori initial condition for the window is the one defined in the
% estimation problem.
x0 = problem.x0;

% The covariance of the intial condition is also the one defined in the
% estimation problem.
Q0 = problem.Q0;

% Let's create now this window using the Window class of Mohotools.
% Note that we also pass the problem to the window, becasue the problem is what
% holds the experimental measurements.
window = mohotools.Window(problem, [k0 k], x0, Q0);

%% 2.b Solve the window.

% Lastly, we can solve the window.
solution = window.solve();

% The solution provides cost function value at the optimum:
solution.J
% , the optimal initial condition:
solution.x0
% , the time vector of the estimation of the states:
solution.t
% , the time evolution of the state vector of the estimation:
solution.x
% , the time evolution of the covariance of the state vector of the estimation:
solution.P 
% , the time evolution of the output:
solution.y

%% 2.c Plot the solution.

% 'utils.plot_data' is a helper function that plots the given dataset.
utils.plot_data('exp01');

% Finally, we can plot the solution of the window to check how it fits the
% experimental data.
subplot(1,2,1);
errorbar(solution.t, solution.x(:,1), 2*sqrt(squeeze(solution.P(:,1,1))), 'o', 'DisplayName', 'S (estimate)');
errorbar(solution.t, solution.x(:,2), 2*sqrt(squeeze(solution.P(:,2,2))), 'o', 'DisplayName', 'X (estimate)');

subplot(1,2,2);
errorbar(solution.t, solution.x(:,1), 2*sqrt(squeeze(solution.P(:,1,1))), 'o-', 'DisplayName', 'S (estimate)');
