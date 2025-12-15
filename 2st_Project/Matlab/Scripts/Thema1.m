% Thema 1 - Plot the function
clear; clc; close all;

f = @(x, y) x.^3 .* exp(-x.^2 -y.^4);

x = linspace(-3, 3, 150);
y = linspace(-3, 3, 150);
[X, Y] = meshgrid(x, y);
F = f(X, Y);

% Plot the function and mark its minimum and maximum points
[min_value, min_index] = min(F(:));
[min_row, min_column] = ind2sub(size(F), min_index);
x_min = X(min_row, min_column);
y_min = Y(min_row, min_column);
[max_value, max_index] = max(F(:));
[max_row, max_column] = ind2sub(size(F), max_index);
x_max = X(max_row, max_column);
y_max = Y(max_row, max_column);

figure;
surf(X, Y, F);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('$f(x, y)$', 'Interpreter', 'latex', 'FontSize', 14);
title('Plot of $f(x, y) = x^3e^{- x^2- y^4}$ and its minimum and maximum', 'Interpreter', 'latex', 'FontSize', 16);
colorbar;
hold on;
plot3(x_min, y_min, min_value, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(x_min, y_min, min_value - 0.1, sprintf('(%.3f, %.3f, %.3f)', x_min, y_min, min_value), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 12);
hold on;
plot3(x_max, y_max, max_value, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(x_max, y_max, max_value + 0.1, sprintf('(%.3f, %.3f, %.3f)', x_max, y_max, max_value), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);
