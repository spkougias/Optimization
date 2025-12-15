% Define the function
f = @(x, y) 1/3*x.^2 + 3*y.^2;
x = linspace(-10, 5, 100);
y = linspace(-8, 12, 100);
[X, Y] = meshgrid(x, y);
F = f(X, Y);

% Plot the function and mark its minimum point
[min_value, min_index] = min(F(:));
[min_row, min_column] = ind2sub(size(F), min_index);
x_min = X(min_row, min_column);
y_min = Y(min_row, min_column);

figure;
surf(X, Y, F);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('$f(x, y)$', 'Interpreter', 'latex', 'FontSize', 14);
title('Plot of $f(x_1, x_2) = \frac{1}{3} x_1^2 + 3x_2^2$ and its minimum point', ...
     'Interpreter', 'latex', 'FontSize', 16);
colorbar;
hold on;
plot3(x_min, y_min, min_value, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
text(x_min, y_min, min_value - 20, sprintf('(%.3f, %.3f, %.3f)', x_min, y_min, min_value), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 12);
