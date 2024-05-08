function [problem] = load_estimation_problem(dataset_name)
% LOAD_ESTIMATION_PROBLEM loads the estimation problem with the given
% dataset.

%% Init an empty estimation problem.

problem = mohotools.EstimationProblem();

%% Load the measurements.

measurements = readtable(['data/' + dataset_name + '/measurements.csv']);
measurements_covariance = readmatrix(['data/' + dataset_name + '/measurements_covariance.csv']);

% Experiment time vector.
problem.t = measurements.time;

% Experiment measurements.
problem.y = measurements.S;

% Covariance matrix of the measurements.
problem.R = measurements_covariance;

%% Load the _a priori_ knowledge of the experiment's initial condition.

% A priori initial condition.
problem.x0 = [12; 0.5];

% Covariance matrix of the a priori initial condition.
problem.Q0 = diag([1 1]);

%% States bounds.

% States upper_bound.
problem.x_ub = [20; 20];

% States lower bound.
problem.x_lb = [0; 0];

%% Load the model.

p = readtable("+model/parameter.csv");

% Function handle which evaluates the model ode.
problem.ode = @(t,x) model.ode(t, x, p.value);

% Function handle which evaluates the model output.
problem.output = @(t,x) model.output(t, x, p.value);

end