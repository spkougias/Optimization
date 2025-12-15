clear; clc; close all;

% Define
syms x y;
f_sym = 1/3*x^2 + 3*y^2;
vars = [x, y];
x0 = [-1, 1];   % starting point
epsilon = 0.001;
max_iters = 1000;

gamma_values = [0.1, 0.3, 3, 5];

% 3. Run Loop
for i = 1:length(gamma_values)
    gamma = gamma_values(i);
    fprintf('Steepest Descent with gamma = %.1f...\n', gamma);
    
    [x_hist, k] = steepest_descent(f_sym, vars, x0, epsilon, gamma, max_iters);

    f_hist = zeros(size(x_hist, 1), 1);
    f_handle = matlabFunction(f_sym, 'Vars', vars);
    for j = 1:size(x_hist, 1)
        f_hist(j) = f_handle(x_hist(j,1), x_hist(j,2));
    end
    
    % Plot convergence of f value
    figure('Name', sprintf('Gamma = %.1f', gamma));
    subplot(2, 1, 1);
    plot(0:k, f_hist, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    title(sprintf('Obj. Function vs Iterations (\\gamma = %.1f)', gamma));
    xlabel('Iteration k');
    ylabel('f(x_k)');
    grid on;
    
    % Plot trajectory of vars
    subplot(2, 1, 2);
    plot(0:k, x_hist(:,1), '-r', 'LineWidth', 1.5); hold on;
    plot(0:k, x_hist(:,2), '-b', 'LineWidth', 1.5);
    legend('x_1', 'x_2');
    title(sprintf('Variables vs Iterations (\\gamma = %.1f)', gamma));
    xlabel('Iteration k');
    ylabel('Value');
    grid on;
    

end
