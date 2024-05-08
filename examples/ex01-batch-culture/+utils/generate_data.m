close all;
clear all;

sigma = [0.2];
exp_ind = 1;

%% Simulate model.

p = readtable("+model/parameter.csv");
x0 = readtable("+model/initial_condition.csv");
x0.value = [10; 1];
ode = @(t,x) model.ode(t, x, p.value);
tspan = 0:1:30;

[t, x] = ode45(ode, tspan, x0.value);

y = [];
for i = 1:length(t)
    y(i,:) = model.output(t(i), x(i,:), p.value)';
end

y_noise = y + sigma.*randn(size(y));

%% Plot simulation.

close all;

figure(1);
hold on;
plot(t, x);
legend(x0.name);
title('Ground truth');

figure(2);
hold on;
plot(t, y);
plot(t, y_noise, 'o');
legend('y', 'measurement');
title('Noisy measurements');

%% Save data.

% Ground truth.
time = t;
S = x(:,1);
X = x(:,2);
T = table(time, S, X);

writetable(T, sprintf("data/exp0%d/ground_truth.csv", exp_ind));

% Measurements.
S = y_noise(:,1);
T = table(time, S);
writetable(T, sprintf("data/exp0%d/measurements.csv", exp_ind));

% Covariance matrix.
covariance = diag(sigma.^2);
writematrix(covariance, sprintf("data/exp0%d/measurements_covariance.csv", exp_ind));