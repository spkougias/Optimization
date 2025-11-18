% Define functions and initial interval
f1 = @(x) 5.^x + (2 - cos(x)).^2;
f2 = @(x) (x - 1).^2 + exp(x - 5).*sin(x + 3);
f3 = @(x) exp(-3*x) - (sin(x - 2) - 2).^2;
a = -1;
b = 3;

% Plot the functions over the given interval and mark their minimum point
interval = [a, b];
x_min = [fminbnd(f1, interval(1), interval(2)), ...
         fminbnd(f2, interval(1), interval(2)), ...
         fminbnd(f3, interval(1), interval(2))];
y_min = [f1(x_min(1)), f2(x_min(2)), f3(x_min(3))];

figure;
subplot(1, 3, 1);
fplot(f1, interval, 'b-', 'LineWidth', 1.5);
hold on;
plot(x_min(1), y_min(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
ylim([0, max(ylim)]);
xlabel('x');
ylabel('f_1(x)');
title('f_1(x)');
grid on;
text(x_min(1), y_min(1) +5, sprintf('(%.3f, %.3f)', x_min(1), y_min(1)), ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

subplot(1, 3, 2);
fplot(f2, interval, 'b-', 'LineWidth', 1.5);
hold on;
plot(x_min(2), y_min(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
ylim([0, max(ylim)]);
xlabel('x');
ylabel('f_2(x)');
title('f_2(x)');
grid on;
text(x_min(2), y_min(2) + 0.8, sprintf('(%.3f, %.3f)', x_min(2), y_min(2)), ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

subplot(1, 3, 3);
fplot(f3, interval, 'b-', 'LineWidth', 1.5);
hold on;
plot(x_min(3), y_min(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
ylim([-5, 100]);
xlabel('x');
ylabel('f_3(x)');
title('f_3(x)');
grid on;
text(x_min(3) - 0.1, y_min(3) + 5, sprintf('(%.3f, %.3f)', x_min(3), y_min(3)), ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

sgtitle('True minimum values of the functions f_1, f_2 and f_3 in [-1, 3]');
