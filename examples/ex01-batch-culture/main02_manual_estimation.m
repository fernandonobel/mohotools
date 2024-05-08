%% Main 02 -- Find the optimal initial condition manually.
%
% The goal of this script is to show:
% 
% 1. How full horizon estimation works in an intuitive way.

clear all;
close all;

%% 1. How full horizon estimation works in an intuitive way.
%
% The fundamental idea of full-horizon estimation is quite simple: find the
% intial condition for the experiment that best explain the avilable
% experimental measurements.
%
% Before using Mohotools methods, let's implement this procedure in a manual
% fashion to ensure that you understand the main idea. Then, in the following
% scritps, we will use Mohotools to improve this manual procedure.

% First let's load the measurements.
measurements = readtable('data/exp01/measurements.csv');

% And plot them.
f = figure();
hold on;
grid on;
plot(measurements.time, measurements.S, 'o', 'DisplayName', 'S (measurement)');
legend('show');
ylabel('Concentration [g/L]');
xlabel('Time [h]');

% Then, let's choose a arbitrary intial condition for the model.
x0 = [12 0.5];

% And simulate and plot the predition of the model for that given initial
% condition.
p_table = readtable("+model/parameter.csv");
p = p_table.value;
ode = @(t,x) model.ode(t, x, p);
tspan = 0:0.1:30;
[t, x] = ode45(ode, tspan, x0);
plot(t, x(:,1), 'DisplayName', 'S (model)');
plot(t, x(:,2), 'DisplayName', 'X (model)');

% The goal is to tune the initial condition so the output of the model (the
% substrate) fits the experimental measurements.
%
% Try to manually modify the value of 'x0' to make that.
