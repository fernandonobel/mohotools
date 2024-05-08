function [] = plot_window(solution, name)

arguments
    solution
    name (1,:) char
end

subplot(1,2,1);
errorbar(solution.t, solution.x(:,1), 2*sqrt(squeeze(solution.P(:,1,1))), 'o-', 'DisplayName', ['S (' name ')']);
errorbar(solution.t, solution.x(:,2), 2*sqrt(squeeze(solution.P(:,2,2))), 'o-', 'DisplayName', ['X (' name ')']);

subplot(1,2,2);
% plot(solution.t, solution.y(:,1), '--', 'DisplayName', 'S (estimate)');
errorbar(solution.t, solution.x(:,1), 2*sqrt(squeeze(solution.P(:,1,1))), 'o-', 'DisplayName', ['S (' name ')']);

end