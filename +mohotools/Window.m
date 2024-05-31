classdef Window < handle
    % WINDOW handles all the information needed to solve an individual
    % window of either full- or moving-horizon estimation algorithms.
    %
    % In particular:
    % - The EstimationProblem which handles the measurements and the model.
    % - The window start and end indices.
    % - The a priori initial condition of the window (and its covariance).
    
    properties
        % The EstimationProblems which handles the experimental data and
        % the model.
        problem (1,1) mohotools.EstimationProblem
        
        % The window interval.
        % It must be a two-element vector [k0 kf] specifing the start and
        % end indices of the window.
        interval (1,2) {mustBeInteger}
        
        % The a priori estimation for the initial condition at the start
        % of the window (k0).
        x0
        
        % The covariance of the a priori estimation of the initial
        % condition (x0).
        Q0
        
        % Time span for the simulation of the solution.
        tspan (1,:) double
    end
    
    properties (Dependent)
        % The time vector inside the window.
        t
        
        % The measurements matrix inside the window.
        y
    end
    
    methods
        function obj = Window(problem, interval, x0, Q0, tspan)
            % Constructor method for Window.
            
            arguments
                problem
                interval
                x0
                Q0
                tspan (1,:) double = [];
            end
            
            obj.problem = problem;
            obj.interval = interval;
            obj.x0 = x0;
            obj.Q0 = Q0;
            
            if isempty(tspan)
                obj.tspan = obj.t;
            else
                obj.tspan = tspan;
            end
        end
        
        function [solution] = solve(obj)
            %% SOLVE finds the optimal initial condition for the window 
            % which best fits the experimental data.
            
            arguments
                obj (1,1) mohotools.Window
            end
            
            %%% Define the nonlinear estimation problem.
            
            % Options for fmincon solver.
            options = optimset();
            options.display = 'iter';
            options.PlotFcns = @optimplotfval;
            options.MaxFunEvals = 1e9;
            
            % Definition of the non-linear estimation problem.
            theta0 = obj.x0;
            ub = obj.problem.x_ub;
            lb = obj.problem.x_lb;
            cost_function = @(theta) obj.cost_function(theta);
            
            %%% Solve the nonlinear estimation problem.
            
            [theta_hat, J] = fmincon(cost_function, theta0, [], [], [], [], lb, ub, [], options);
            
            %%% Simulate the solution.
            
            % Get initial condition.
            x0_hat = theta_hat;
                       
            % Simulate the model.
            if length(obj.tspan) >= 2
                [t, x] = ode45(@(t,x) obj.problem.ode(t,x), obj.tspan, x0_hat);
            else
                t = obj.tspan;
                x = x0_hat';
            end
            
            if length(obj.tspan) == 2
                t = [t(1); t(end)];
                x = [x(1,:); x(end,:)];
            end
            
            % Calculate the output.
            y = [];
            for i = 1:length(t)
                y(i,:) = obj.problem.output(t(i), x(i,:));
            end
            
            % Calculate the conficence interval of the solution.
            [P] = obj.calculate_confidence_interval(x0_hat);
            
            solution.x0 = theta_hat;
            solution.J = J;
            solution.t = t;
            solution.x = x;
            solution.y = y;
            solution.P = P;
        end
        
        function [J] = cost_function(obj, theta)
            %% COST_FUNCTION
            
            arguments
                % The window object.
                obj (1,1) mohotools.Window
                % The decision vector.
                theta (:,1) double
            end
            
            % Initialize cost value.
            J = 0;
            
            % Initial condition for the window.
            x0 = theta;
            
            % Simulation time span.
            tspan = obj.t;
            
            % Simulate the model.
            if length(tspan) >= 2
              [t, x] = ode45(@(t,x) obj.problem.ode(t,x), tspan, x0);
            else
               t = tspan;
               x = x0';
            end
            
            % 
            if length(tspan) == 2
                t = [t(1); t(end)];
                x = [x(1,:); x(end,:)];
            end
            
            % Add the arrival cost.
            e = obj.x0 - x0;
            J = J + e' * inv(obj.Q0) * e;

            % Evaluate the fitting of the model prediciton to the experimental data.
            for i = 1:length(t)
                y_model = obj.problem.output(t(i), x(i,:));
                e = obj.y(i,:)' - y_model;
                
                not_nan = ~isnan(e);
                e = e(not_nan);
                
                J = J + e' * inv(obj.problem.R(not_nan, not_nan)) * e;
            end
            
        end
        
        function [P] = calculate_confidence_interval(obj, x0)
            %% CALCULATE_CONFIDENCE_INTERVAL This functions calculates the confidence
            % interval of estimation of the states for a given solution of a window.
            
            arguments
                % The window object.
                obj (1,1) mohotools.Window
                % Estimated initial condition for the window.
                x0 (:,1) double
            end
            
            % Number of states.
            N_states = length(x0);
            
            %%% Define symbolic variables and equations.
            
            % Time.
            t = sym('t');
            
            % State vector.
            x = sym('x', [N_states 1]);
            
            % The ode function.
            f = obj.problem.ode(t, x);
            
            % The ouput function.
            h = obj.problem.output(t, x);
            
            % The jacobian of the ode function with respect the states.
            df_dx = jacobian(f, x);
            
            % The jacobian of the output function with respect the states.
            dh_dx = jacobian(h, x);
            dh_dx_function = matlabFunction(dh_dx, "Vars", {t, x});
            
            %%% Calculate the sensitivity functions.
            
            % Sensitivity functions.
            G = sym('G', [N_states N_states]);
            
            % The ode function for the sensitivity functions.
            dG_dt = df_dx*G;
            
            % Initial condition of sensitivity functions.
            G0 = eye(N_states, N_states);
            
            % Create extended system (states + sensitivities).
            z = [
                x
                reshape(G, N_states*N_states, 1)
                ];
            
            z0 = [
                x0
                reshape(G0, N_states*N_states, 1)
                ];
            
            dz_dt = [
                f
                reshape(dG_dt, N_states*N_states, 1)
                ];
            
            dz_dt_function = matlabFunction(dz_dt, "Vars", {t, z});
            
            % Simulate the extended model.
            if length(obj.t) >= 2
                [t, z] = ode45(dz_dt_function, obj.t, z0);
            else
                t = obj.t;
                z = z0';
            end
            
            % Get the states.
            x = z(:, 1:N_states);
            
            % Get the sensitivities.
            G = z(:, N_states+1:end);
            
            %%% Calculate the covariance of the initial condition.
            
            P0 = zeros(N_states, N_states);
            
            % Add the information of the arrival cost.
            P0 = inv(obj.Q0);
            
            % Add the information from the measurements.
            for k = 1:length(obj.t)
                C_k = dh_dx_function(t(k), x(k,:));
                G_k = reshape(G(k,:), N_states, N_states);
                
                P0 = P0 + (C_k*G_k)' * inv(obj.problem.R) * (C_k*G_k);
            end
            
            P0 = inv(P0);
            
            %%% Calculate the covariance of the evolution of the state.
            
            % Simulate the extended model.
            if length(obj.tspan) >= 2
                [t, z] = ode45(dz_dt_function, obj.tspan, z0);
            else
                t = obj.tspan;
                z = z0';
            end
            
            % Get the states.
            x = z(:, 1:N_states);
            
            % Get the sensitivities.
            G = z(:, N_states+1:end);
            
            P = zeros(length(obj.tspan), N_states, N_states);
            
            for k = 1:length(obj.tspan)
                G_k = reshape(G(k,:), N_states, N_states);
                P_k = G_k * P0 * G_k';
                P(k,:,:) = P_k;
            end
            
        end
        
        function value = get.t(obj)
            %% Get the time vector of the window.
            
            % Get the indices of measurements inside the window.
            ind = obj.interval(1) : obj.interval(2);
            value = obj.problem.t(ind);
        end
        
        function value = get.y(obj)
            %% Get measurements inside the window.
            
            % Get the indices of measurements inside the window.
            ind = obj.interval(1) : obj.interval(2);
            value = obj.problem.y(ind, :);
        end
    end
end
