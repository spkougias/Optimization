% Thema 4
clear; clc; close all;

syms x y;
f_sym = x^3 * exp(-x^2 - y^4);
vars = [x, y];
f_handle = matlabFunction(f_sym);

epsilon = 0.001;
start_points = [0, 0; -1, -1; 1, 1];

strategies = {
    'Fixed Step (gamma=0.1)', 1, 0.1; 
    'Exact Search', 3, 2.5; 
    'Armijo Rule', 2, 1.0
};

for i = 1:size(start_points, 1)
    x0 = start_points(i, :);
    
    for s = 1:size(strategies, 1)
        strat_name = strategies{s, 1};
        mode = strategies{s, 2};
        gamma_val = strategies{s, 3};
        
        fprintf('Start (%.1f, %.1f) - %s... \n', x0(1), x0(2), strat_name);
        
        try
            [x_hist, iters] = levenberg_marquardt(f_sym, vars, x0, epsilon, gamma_val, mode);
        catch
            x_hist = [x0; x0]; iters = 0;
        end
        
        f_vals = zeros(size(x_hist, 1), 1);
        for k = 1:size(x_hist, 1)
            f_vals(k) = f_handle(x_hist(k,1), x_hist(k,2));
        end
        
        % Plotting
        fig_name = sprintf('LM - Start (%.0f, %.0f) - %s', x0(1), x0(2), strat_name);
        figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [100, 100, 1000, 400]);
        
        subplot(1, 2, 1); 
        plot(0:length(f_vals)-1, f_vals, '-or', 'LineWidth', 1.5, 'MarkerSize', 4);
        title(['Convergence: ' strat_name]); xlabel('k'); ylabel('f(x)'); grid on;
        
        subplot(1, 2, 2);
        margin = 0.5;
        x_min = min([-2.5; x_hist(:,1)]) - margin;
        x_max = max([ 2.5; x_hist(:,1)]) + margin;
        y_min = min([-2.5; x_hist(:,2)]) - margin;
        y_max = max([ 2.5; x_hist(:,2)]) + margin;
        
        xv = linspace(x_min, x_max, 100);
        yv = linspace(y_min, y_max, 100);
        [X, Y] = meshgrid(xv, yv);
        Z = f_handle(X, Y);
        
        contour(X, Y, Z, 30, 'LineColor', [0.7 0.7 0.7]); hold on;
        plot(x_hist(:,1), x_hist(:,2), 'b.-', 'LineWidth', 1.2, 'MarkerSize', 8);
        plot(x0(1), x0(2), 'gx', 'LineWidth', 2, 'MarkerSize', 10);
        plot(x_hist(end,1), x_hist(end,2), 'ro', 'LineWidth', 2, 'MarkerSize', 10);
        
        title(['Path: ' strat_name]);
        axis([x_min x_max y_min y_max]); grid on;
    end
end