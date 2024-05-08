function [dx_dt] = ode(t, x, p)
% ODE evaluates the derivative vector.

arguments
    % Current time in the simulation.
    t (1,1)
    % State column-vector.
    x (:,1)
    % Parameter column-vector.
    p (:,1)
end

%% Define the states.

S = x(1);
X = x(2);

%% Define the parameters.

mu_max = p(1);
K = p(2);

%% Evaluate the equations of the model.

phi = mu_max * S/(K+S) * X;

dS_dt = -phi;
dX_dt = phi;

%% Return derivative vector.

dx_dt = [
    dS_dt
    dX_dt
    ];

end
