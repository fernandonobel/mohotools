%% Main 00 -- Problem definition
%
% The goal of this script is to define the problem.

%% Monod kinetics in a batch culture
% 
% Let's consider a batch bioreactor where we only have the following reaction:
% 
%  [S] -> [X] ; phi(S,X)
% 
% where:
% 
% - 'S' is the substrate concentration [g/L]
% - 'X' is the cell concentration [g/L]
% - 'phi(S,X)' is the substrate consumption reaction rate, defined as:
% 
%   phi(S,X) = mu_max * S/(K+S) * X
% 
% with 'mu_max=0.2' and 'K=1'.
% 
% Applying mass action kinetics, we can obtain the following ODE system:
% 
%   der(S) = - phi(S,X)
%   der(X) = phi(S,X) 
% 
% - Measurements of 'S' are available every half hour with a zero-mean Gaussian 
%   noise with stantard deviation of $0.2$. 
% - The initial condition of the experiment is 'S_0=10' and 'X=1' (but we assume
%   that is unknown for the observer).
% - The a priori estimation of the initial condition is 'S_0=12' and 'X=0.5'.
% - The model parameters we have identified are: 'mu_max=0.1958' and 'K=0.8929'.
% 
% Design a moving horizon estimator to estimate the biomass using the substrate
% measurements.

%% References
%
% This example is from this slides: 
%
%   "Introduction to linear and nonlinear parameter identification - with 
%   Matlab worked examples. Alain Vande Wouwer. 2023."
