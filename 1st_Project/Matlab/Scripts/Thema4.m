clc
% --- Thema 4 ---

% --- Setup ---
syms x;

f1 = @(x) 5^x + (2 - cos(x))^2;
f2 = @(x) (x - 1)^2 + exp(x - 5)*sin(x + 3);
f3 = @(x) exp(-3*x) - (sin(x - 2) - 2)^2;
a = -1;
b = 3;

% Find the derivative of the functions
df1_sym = diff(f1, x);
df2_sym = diff(f2, x);
df3_sym = diff(f3, x);
df1 = matlabFunction(df1_sym);
df2 = matlabFunction(df2_sym);
df3 = matlabFunction(df3_sym);

clear x;

% --- Task 1: n vs l ---
l_values = linspace(0.001, 0.5, 30);
n_values = zeros(3,length(l_values));
for l_idx = 1:length(l_values)
    [~,n_values(1,l_idx)] = bisection_method_with_derv(df1, a, b, l_values(l_idx));
    [~,n_values(2,l_idx)] = bisection_method_with_derv(df2, a, b, l_values(l_idx));
    [~,n_values(3,l_idx)] = bisection_method_with_derv(df3, a, b, l_values(l_idx));
end
figure;
sgtitle('Bisection (Derivative): Number of Evaluations vs. l (accuracy)');
subplot(1, 3, 1);
plot(l_values, n_values(1,:), 'b-o', 'LineWidth', 1.5);
title('f_1(x)'); xlabel('l (accuracy)'); ylabel('n (evaluations)'); grid on;
subplot(1, 3, 2);
plot(l_values, n_values(2,:), 'r-o', 'LineWidth', 1.5);
title('f_2(x)'); xlabel('l (accuracy)'); ylabel('n (evaluations)'); grid on;
subplot(1, 3, 3);
plot(l_values, n_values(3,:), 'g-o', 'LineWidth', 1.5);
title('f_3(x)'); xlabel('l (accuracy)'); ylabel('n (evaluations)'); grid on;

% --- Task 2: [ak, bk] vs k (various l) ---
l_vector = [0.005, 0.01, 0.1, 0.5];
intervals_f1 = cell(1, length(l_vector));
intervals_f2 = cell(1, length(l_vector));
intervals_f3 = cell(1, length(l_vector));
for l_idx = 1:length(l_vector)
    [intervals_f1{l_idx}, ~] = bisection_method_with_derv(df1, a, b, l_vector(l_idx));
    [intervals_f2{l_idx}, ~] = bisection_method_with_derv(df2, a, b, l_vector(l_idx));
    [intervals_f3{l_idx}, ~] = bisection_method_with_derv(df3, a, b, l_vector(l_idx));
end
figure;
sgtitle('Bisection (Derivative): Interval Endpoints vs. Iteration k');
colors = lines(length(l_vector));
legends = cell(1, length(l_vector));
for l_idx = 1:length(l_vector)
    legends{l_idx} = ['l = ', num2str(l_vector(l_idx), '%.3f')];
end
% F1
subplot(1, 3, 1);
hold on;
for l_idx = 1:length(l_vector)
    k_values = 1:size(intervals_f1{l_idx}, 1);
    plot(k_values, intervals_f1{l_idx}(:,1), '-o', 'Color', colors(l_idx, :), 'LineWidth', 1.3, 'MarkerSize', 6, 'DisplayName', legends{l_idx});
    plot(k_values, intervals_f1{l_idx}(:,2), '-x', 'Color', colors(l_idx, :), 'LineWidth', 1.3, 'MarkerSize', 6, 'HandleVisibility', 'off');
end
title('f_1(x)'); xlabel('k (iteration)'); ylabel('[a_k, b_k]');
legend('show', 'Location', 'best'); grid on; hold off;
% F2
subplot(1, 3, 2);
hold on;
for l_idx = 1:length(l_vector)
    k_values = 1:size(intervals_f2{l_idx}, 1);
    plot(k_values, intervals_f2{l_idx}(:,1), '-o', 'Color', colors(l_idx, :), 'LineWidth', 1.3, 'MarkerSize', 6, 'DisplayName', legends{l_idx});
    plot(k_values, intervals_f2{l_idx}(:,2), '-x', 'Color', colors(l_idx, :), 'LineWidth', 1.3, 'MarkerSize', 6, 'HandleVisibility', 'off');
end
title('f_2(x)'); xlabel('k (iteration)'); ylabel('[a_k, b_k]');
legend('show', 'Location', 'best'); grid on; hold off;
% F3
subplot(1, 3, 3);
hold on;
for l_idx = 1:length(l_vector)
    k_values = 1:size(intervals_f3{l_idx}, 1); 
    plot(k_values, intervals_f3{l_idx}(:,1), '-o', 'Color', colors(l_idx, :), 'LineWidth', 1.3, 'MarkerSize', 6, 'DisplayName', legends{l_idx});
    plot(k_values, intervals_f3{l_idx}(:,2), '-x', 'Color', colors(l_idx, :), 'LineWidth', 1.3, 'MarkerSize', 6, 'HandleVisibility', 'off');
end
title('f_3(x)'); xlabel('k (iteration)'); ylabel('[a_k, b_k]');
legend('show', 'Location', 'best'); grid on; hold off;