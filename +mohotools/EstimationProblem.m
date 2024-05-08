classdef EstimationProblem < handle
    % ESTIMATION_PROBLEM stores all the information of a given estimation
    % problem.
    %
    % In particular:
    % - The experimental data (and its covariance).
    % - The a priori initial condition of the experiment (and its covariance).
    % - The upper and lower bounds of the state vector.
    % - The model ode and output functions.
    
    properties
        %% Experimental measurements.
        
        % Time vector of the measurements.
        t (:,1) double
        
        % Measurements value matrix.
        y (:,:) double
        
        % Covariance matrix of the measurements.
        R (:,:) double
        
        %% A priori initial condition of the window.
        
        % A priori estimation of the initial condition of the experiment.
        x0 (:,1) double
        
        % Covariance matrix of x_0.
        Q0 (:,:) double
        
        %% Upper and lower bounds of the state vector.
        
        % Upper-bound for the state vector.
        x_ub (:,1) double
        
        % Lower-bound for the state vector.
        x_lb (:,1) double
        
        %% Model of the process.
        
        % A funtion handle that evaluates the model ODE:
        %
        % The arguments must be:
        % - t (1)   Simulation time.
        % - x (:,1) State vector.
        % And it must return a column vector with the derivatives.
        ode function_handle = @(t,x) []
        
        % Function that evaluates the output of the model:
        %
        % The arguments must be:
        % - t (1)   Simulation time.
        % - x (:,1) State vector.
        % And it must return a column vector with the output.
        output function_handle = @(t,x) []   
    end
end
