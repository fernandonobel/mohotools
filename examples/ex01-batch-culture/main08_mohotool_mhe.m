%% Main 08 -- Mohotool's method for Moving Horizon Estimation.
%
% The goal of this script is to show:
%
% 1. How can we use Mohotool to easily implement moving-horizon estimation.

close all;
clear all;

%% 1. How can we use Mohotool to easily implement moving-horizon estimation.
%
% In previous scripts, we used the Window class of Mohotools to implement full-
% and moving-horizon estimation. This was done just for didactically reasons. 
% The expected way of using Mohotools is not by manually manipulating each
% window. The normal use of this toolbox is:
%
% 1. The user defines the estimation problem.
% 2. The user applies one of the built-in methods to solve the problem.
%
% In this script, we show how to do moving-horizon estimation this way.

%% Load the estimation problem.

data = "exp01";
problem = utils.load_estimation_problem(data);

%% Solve the problem with Moving Horizon Estimation.

N = 5;
[t, x, P] = mohotools.moving_horizon_estimation(problem, N);

%% Plot result.

figure();

utils.plot_data(data);

subplot(1,2,1);
errorbar(t, x(:,1), 2*sqrt(squeeze(P(:,1,1))), 'o-', 'DisplayName', 'S (mhe)');
errorbar(t, x(:,2), 2*sqrt(squeeze(P(:,2,2))), 'o-', 'DisplayName', 'X (mhe)');

subplot(1,2,2);
errorbar(t, x(:,1), 2*sqrt(squeeze(P(:,1,1))), 'o-', 'DisplayName', 'S (mhe)');
