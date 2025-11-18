clc
% --- Thema 1 ---

% --- Setup ---
f1 = @(x) 5^x + (2 - cos(x))^2;
f2 = @(x) (x - 1)^2 + exp(x - 5)*sin(x + 3);
f3 = @(x) exp(-3*x) - (sin(x - 2) - 2)^2;
a = -1;
b = 3;


% --- Task 1 n vs epsilon (fixed l) --- 
l_fixed = 0.01;

e_values = linspace(1e-5, l_fixed/2 - 1e-5, 30);
n_values_e = zeros(3, length(e_values));

for e = 1:length(e_values)
    [~, n_values_e(1,e)] = bisection_method(f1, a, b, e_values(e), l_fixed);
    [~, n_values_e(2,e)] = bisection_method(f2, a, b, e_values(e), l_fixed);
    [~, n_values_e(3,e)] = bisection_method(f3, a, b, e_values(e), l_fixed);
end

figure;
sgtitle(['Number of Evaluations vs. \epsilon for l = ', num2str(l_fixed, '%.2f')]);

subplot(1, 3, 1);
plot(e_values, n_values_e(1,:), 'b-o', 'LineWidth', 1.5);
title('f_1(x)'); xlabel('\epsilon'); ylabel('n (evaluations)'); grid on;

subplot(1, 3, 2);
plot(e_values, n_values_e(2,:), 'r-o', 'LineWidth', 1.5);
title('f_2(x)'); xlabel('\epsilon'); ylabel('n (evaluations)'); grid on;

subplot(1, 3, 3);
plot(e_values, n_values_e(3,:), 'g-o', 'LineWidth', 1.5);
title('f_3(x)'); xlabel('\epsilon'); ylabel('n (evaluations)'); grid on;

% --- Task 2 n vs l (fixed epsilon) --- 
e_fixed = 0.001;

l_values = linspace(2*e_fixed + 1e-5, 0.5, 30);
n_values_l = zeros(3, length(l_values));

for l_idx = 1:length(l_values)
    [~, n_values_l(1,l_idx)] = bisection_method(f1, a, b, e_fixed, l_values(l_idx));
    [~, n_values_l(2,l_idx)] = bisection_method(f2, a, b, e_fixed, l_values(l_idx));
    [~, n_values_l(3,l_idx)] = bisection_method(f3, a, b, e_fixed, l_values(l_idx));
end

figure;
sgtitle(['Number of Evaluations vs. l for \epsilon = ', num2str(e_fixed, '%.3f')]);

subplot(1, 3, 1);
plot(l_values, n_values_l(1,:), 'b-o', 'LineWidth', 1.5);
title('f_1(x)'); xlabel('l (accuracy)'); ylabel('n (evaluations)'); grid on;

subplot(1, 3, 2);
plot(l_values, n_values_l(2,:), 'r-o', 'LineWidth', 1.5);
title('f_2(x)'); xlabel('l (accuracy)'); ylabel('n (evaluations)'); grid on;

subplot(1, 3, 3);
plot(l_values, n_values_l(3,:), 'g-o', 'LineWidth', 1.5);
title('f_3(x)'); xlabel('l (accuracy)'); ylabel('n (evaluations)'); grid on;

% --- Task 3 [ak, bk] vs k (l varies) --- 
l_vector = [0.005, 0.01, 0.1, 0.5];
intervals_f1 = cell(1, length(l_vector));
intervals_f2 = cell(1, length(l_vector));
intervals_f3 = cell(1, length(l_vector));

for l_idx = 1:length(l_vector)
    [intervals_f1{l_idx}, ~] = bisection_method(f1, a, b, e_fixed, l_vector(l_idx));
    [intervals_f2{l_idx}, ~] = bisection_method(f2, a, b, e_fixed, l_vector(l_idx));
    [intervals_f3{l_idx}, ~] = bisection_method(f3, a, b, e_fixed, l_vector(l_idx));
end

figure;
sgtitle(['Interval Endpoints vs. Iteration k for \epsilon = ', num2str(e_fixed, '%.3f')]);
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