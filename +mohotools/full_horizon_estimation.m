function [t, x, P] = full_horizon_estimation(problem)
%% FULL_HORIZON_ESTIMATION solves a problem using Full Horizon Estimation.

arguments
    problem (1,1) mohotools.EstimationProblem
end


N_measurements = length(problem.t);
N_states = length(problem.x0);

t = problem.t;
x = zeros(N_measurements, N_states);
P = zeros(N_measurements, N_states, N_states);

for k = 1 : N_measurements
    
    % In FHE, 'k0' is always the first measurement.
    k0 = 1;
    
    % A priori initial condition of the window 'k'.
    x_k0 = problem.x0;
    
    % Covariance of the a priori initial condition of the window 'k'.
    Q_k0 = problem.Q0;
    
    % Definition of window 'k'.
    window = mohotools.Window(problem, [k0 k], x_k0, Q_k0);
    
    % Solve window 'k'.
    solution = window.solve();
    
    x(k,:) = solution.x(k, :);
    P(k,:,:) = solution.P(k, :, :);
    
end


end

