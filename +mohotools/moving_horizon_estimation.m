function [t, x, P] = moving_horizon_estimation(problem, N)
%% FULL_HORIZON_ESTIMATION solves a problem using Full Horizon Estimation.

arguments
    problem (1,1) mohotools.EstimationProblem
    N (1,1) {mustBeInteger}
end


N_measurements = length(problem.t);
N_states = length(problem.x0);

t = problem.t;
x = zeros(N_measurements, N_states);
P = zeros(N_measurements, N_states, N_states);

% Cell array of the solutions of each window.
solution = {};

%% Start-up with Full Horizon Estimation (FHE).

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
    solution{k} = window.solve();
    
    x(k,:) = solution{k}.x(k, :);
    P(k,:,:) = solution{k}.P(k, :, :);
end

%% Moving Horizon Estimation (MHE).

for k = N+2 : length(problem.t)
    
    % In MHE, 'k0' is N measurements behind 'k'. 
    k0 = k-N;
    
    % A priori initial condition of the window 'k'.
    x_k0 = solution{k-N-1}.x(end,:)';
    
    % Covariance of the a priori initial condition of the window 'k'.
    Q_k0 = squeeze(solution{k-N-1}.P(end,:,:));
    
    % Time vector for simulating the solution.
    tspan = [k0:k+1] - 1;
    
    % Definition of window 'k'.
    window = mohotools.Window(problem, [k0 k], x_k0, Q_k0, tspan);
    
    % Solve window 'k'.
    solution{k} = window.solve();
    
    x(k,:) = solution{k}.x(end-1, :);
    P(k,:,:) = solution{k}.P(end-1, :, :);
end


end

