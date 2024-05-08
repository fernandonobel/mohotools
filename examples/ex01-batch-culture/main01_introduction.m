%% Main 01 -- Introduction to the standards for the data and the model.
%
% The goal of this script is to show:
%
% 1. How to load the data from a given experiment.
% 2. How to simulate the model for a given initial condition.

clear all;
close all;

%% 1. How to load the data from a given experiment.
%
% In this example, we have data from three different experiments: 'exp01',
% 'exp02' and 'exp03'. These datasets are saved in the 'data' folder and each
% one is saven in its corresponding folder.
%
% The data can be stored using different data standards. For example, the
% data could have been saved as a '.mat' file or, even, as a Excel sheet.
% In this case, the data format used in this case is CSV (comma-separated 
% value).
%
% For each experiment we have three data files:
%
% 1. 'ground_truth.csv'
%
%    This file stores the real time-evolution of the states of the process
%    without noise. This data is normally not available in a real experiment.
%    However, we have provided this data to facilitate undertanding and
%    validating how full- and moving-horizon estimation works.
%
% 2. 'measurements.csv'
%
%    This file stores the measurements performed during the measurements and
%    the time vector of the measurements.
%
% 3. 'measurements_covariance.csv'
%
%    This file stores the covariance matrix of the measurements. In this case,
%    we assume that the measurements are done syncronously.
%
% MATLAB provides functions to read and import CSV files. Then, we can import
% the data of 'exp01' as follows:

states = readtable('data/exp01/ground_truth.csv');
measurements = readtable('data/exp01/measurements.csv');
covariance = readmatrix('data/exp01/measurements_covariance.csv');

% Once we import the dataset, we can use normal MATLAB funtions to manipulate
% and plot the data.

% Plot the ground truth and measurements data.
f = figure();

hold on;
grid on;

plot(states.time, states.S, 'DisplayName', 'S (ground truth)');
plot(states.time, states.X, 'DisplayName', 'X (ground truth)');
plot(measurements.time, measurements.S, 'o', 'DisplayName', 'S (measurement)');

legend('show');
ylabel('Concentration [g/L]');
xlabel('Time [h]');

% Show the value of the covariance matrix.
disp(covariance);

%% 2. How to simulate the model for a given initial condition.
%
% The last part we need is to define the mathematical model we will use for
% full- and moving-horizon estimation. Similar to the experiments, we can use
% different data standards to define the equations and the information of the
% model.
%
% A state-space model is based on two functions: the ODE and the output
% functions. The ODE function calculates the vector of the state derivative and the
% output function calculates the output vector. Both functions takes as an
% input the given time value, the state vector value at that time, and the
% paramer vector value.
% 
% These functions are defined as MATLAB funcitions inside the '+model' folder.
%
% Lastly, the default parameter vector is also saved in '+model' as a CSV files.

% Load the initial condition vector.
x0 = [12 0.5];

% Load the parameter vector.
p_table = readtable("+model/parameter.csv");
p = p_table.value;

% Prepare the function handle needed for the 'ode45' solver.
ode = @(t,x) model.ode(t, x, p);

% The simulation time vector.
tspan = 0:0.1:30;

% Simulate the model.
[t, x] = ode45(ode, tspan, x0);

% Plot simulation result.
plot(t, x(:,1), 'DisplayName', 'S (model)');
plot(t, x(:,2), 'DisplayName', 'X (model)');
