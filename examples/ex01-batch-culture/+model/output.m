function [y] = output(t, x, p)
% OUTPUT evaluates the output vector.

arguments
    % Current time in the simulation.
    t (1,1)
    % State column-vector.
    x (:,1)
    % Parameter column-vector.
    p (:,1)
end

%% Definition of the states.

S = x(1);
X = x(2);

%% Definition of the parameters.

mu_max = p(1);
K = p(2);

%% Evaluate the output vector.

y = [
    S
    ];

end
