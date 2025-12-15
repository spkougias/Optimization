clear; clc; close all;
% Initialize
syms x y;
f = 1/3*x^2 + 3*y^2;
vars = [x, y];
epsilon = 0.01;
a = [-10, -8];
b = [5, 12];

% Thema 2
fprintf('Thema 2...\n');
x0_p2 = [5, -5];
s_p2 = 5;
gamma_p2 = 0.5;
[x_hist2, k2] = projected_steepest_descent(f, vars, x0_p2, epsilon, gamma_p2, s_p2, a, b, false);
plot_results(x_hist2, f, vars, 2, 'Thema 2');

% Thema 3
fprintf('Thema 3...\n');
x0_p3 = [-5, 10];
s_p3 = 15;
gamma_p3 = 0.1;
% Note: We set diminish_s = true here to fix the oscillation
[x_hist3, k3] = projected_steepest_descent(f, vars, x0_p3, epsilon, gamma_p3, s_p3, a, b, true);
plot_results(x_hist3, f, vars, 3, 'Thema 3');

% Thema 4
fprintf('Thema 4...\n');
x0_p4 = [8, -10]; % Infeasible start
s_p4 = 0.1;
gamma_p4 = 0.2;
[x_hist4, k4] = projected_steepest_descent(f, vars, x0_p4, epsilon, gamma_p4, s_p4, a, b, false);
plot_results(x_hist4, f, vars, 4, 'Thema 4');

% --- Local Plotting Function ---
function plot_results(x_vector, f_sym, vars, fig_num, title_str)
    % Helper to plot f values and x coordinates
    f_vals = zeros(size(x_vector, 1), 1);
    for i = 1:size(x_vector, 1)
        f_vals(i) = double(subs(f_sym, vars, x_vector(i, :)));
    end

    figure(fig_num);
    subplot(2,1,1);
    plot(1:size(x_vector,1), f_vals, '-o', 'LineWidth', 1.5);
    title([title_str ' - Objective Function']);
    xlabel('k'); ylabel('f(x)'); grid on;

    subplot(2,1,2);
    plot(1:size(x_vector,1), x_vector(:,1), '-o', 'LineWidth', 1.5); hold on;
    plot(1:size(x_vector,1), x_vector(:,2), '-x', 'LineWidth', 1.5);
    legend('x_1', 'x_2');
    title([title_str ' - Coordinates']);
    xlabel('k'); ylabel('Value'); grid on;
end