function [] = plot_data(name)

arguments
    name (1,:) char
end

ground_truth = readtable(['data/' name '/ground_truth.csv']);
measurements = readtable(['data/' name '/measurements.csv']);
covariance = readmatrix(['data/' name '/measurements_covariance.csv']);

sigma = sqrt(diag(covariance)).*ones(size(measurements.time));

sgtitle(['Data: "' name '"']);

%% Plot ground truth.

subplot(1,2,1);

hold on;
grid on;

plot(ground_truth.time, ground_truth.S, 'DisplayName', 'S (ground truth)');
plot(ground_truth.time, ground_truth.X, 'DisplayName', 'X (ground truth)');

legend('show');
title('States');
ylabel('Concentration [g/L]');
xlabel('Time [h]');

%% Plot measurements.

subplot(1,2,2);

hold on;
grid on;

errorbar(measurements.time, measurements.S, 2*sigma(:,1), 'o', 'DisplayName', 'S (measurement)');

legend('show');
title('Output');
ylabel('Concentration [g/L]');
xlabel('Time [h]');

end