% work1_methods.m
% Clear workspace and set up functions and parameters for the assignment.
clear; close all; clc;

%% === 1. Define the functions to minimize (edit if needed) ===
% NOTE: If the exact mathematical expressions in the PDF differ, replace these.
f{1} = @(x) 5 + (2 - cos(x)).^2;             % f1(x)
f{2} = @(x) (x - 1).^2 + exp(1).*sin(x + 3); % f2(x)
f{3} = @(x) exp( - (sin(x - 2) - 2) );       % f3(x)

% Optional analytic derivatives (if you want to use derivative-based method
% supply df{1}, df{2}, df{3}. If you prefer numerical derivatives, set to [].
df{1} = @(x) 2*(2 - cos(x)).*sin(x);                % derivative of f1
df{2} = @(x) 2*(x - 1) + exp(1).*cos(x + 3);        % derivative of f2
df{3} = @(x) exp( - (sin(x - 2) - 2) ) .* ( - cos(x - 2) ); % derivative of f3

%% === 2. Common parameters ===
a0 = -1; b0 = 3;             % initial interval from the assignment
default_eps = 1e-3;         % epsilon for the dichotomy method
default_l = 0.01;           % example final interval length

%% === 3. Implementations (we call these below). 
% Each returns: a,b,iteration_count,f_eval_count,history (history stores [a b] per iter)

% --- Experiment 1 (Dichotomy) : vary eps with fixed l=0.01 ---
eps_list = [1e-1, 1e-2, 1e-3, 1e-4];
l_fixed = 0.01;
results_dich_eps = cell(3,1);
figure('Name','Dichotomy: function evaluations vs epsilon');
hold on;
for i=1:3
    fcounts = zeros(size(eps_list));
    for j=1:length(eps_list)
        [af,bf,it,fcount,~] = dichotomy(f{i}, a0, b0, l_fixed, eps_list(j));
        fcounts(j) = fcount;
    end
    plot(log10(eps_list), fcounts, '-o', 'DisplayName', sprintf('f_%d',i));
    results_dich_eps{i} = fcounts;
end
xlabel('log10(epsilon)'); ylabel('Number of f evaluations');
title('Dichotomy: f-evals vs epsilon (l = 0.01)');
legend('Location','best'); grid on; hold off;

% --- Experiment 2 (Dichotomy) : vary l with fixed eps=1e-3 ---
l_list = [1e-1, 5e-2, 1e-2, 5e-3, 1e-3];
eps_fixed = 1e-3;
figure('Name','Dichotomy: function evaluations vs l');
hold on;
for i=1:3
    fcounts = zeros(size(l_list));
    for j=1:length(l_list)
        [af,bf,it,fcount,~] = dichotomy(f{i}, a0, b0, l_list(j), eps_fixed);
        fcounts(j) = fcount;
    end
    plot(log10(l_list), fcounts, '-o', 'DisplayName', sprintf('f_%d',i));
    results_dich_eps{i} = fcounts;
end
xlabel('log10(l)'); ylabel('Number of f evaluations');
title(sprintf('Dichotomy: f-evals vs l (eps = %.1e)', eps_fixed));
legend('Location','best'); grid on; hold off;

% --- Endpoint trace (a_k and b_k) for various l values (Dichotomy) ---
l_trace = [0.1, 0.05, 0.01];
for i=1:3 % for each function
    figure('Name',sprintf('Dichotomy endpoints f_%d',i));
    for j=1:length(l_trace)
        [a_fin,b_fin,iter,fcount,history] = dichotomy(f{i}, a0, b0, l_trace(j), eps_fixed);
        % history is N x 2 matrix (a_k, b_k)
        subplot(length(l_trace),1,j);
        k = 0:size(history,1)-1;
        plot(k, history(:,1), '-o', k, history(:,2), '-x');
        xlabel('iteration k'); ylabel('interval endpoints');
        title(sprintf('f_%d, l = %.3g, iterations=%d', i, l_trace(j), iter));
        legend('a_k','b_k','Location','best');
        grid on;
    end
end

%% === 4. Golden Section tests: vary l and plot f-evals ===
l_list_g = [1e-1, 5e-2, 1e-2, 5e-3, 1e-3];
figure('Name','Golden section: f-evals vs l');
hold on;
for i=1:3
    fcounts = zeros(size(l_list_g));
    for j=1:length(l_list_g)
        [af,bf,it,fcount,~] = golden(f{i}, a0, b0, l_list_g(j));
        fcounts(j) = fcount;
    end
    plot(log10(l_list_g), fcounts, '-o', 'DisplayName', sprintf('f_%d',i));
end
xlabel('log10(l)'); ylabel('Number of f evaluations');
title('Golden section: f-evals vs l');
legend('Location','best'); grid on; hold off;

%% === 5. Fibonacci method: tests (vary l) ===
figure('Name','Fibonacci: f-evals vs l');
hold on;
for i=1:3
    fcounts = zeros(size(l_list_g));
    for j=1:length(l_list_g)
        [af,bf,it,fcount,~] = fibonacci(f{i}, a0, b0, l_list_g(j));
        fcounts(j) = fcount;
    end
    plot(log10(l_list_g), fcounts, '-o', 'DisplayName', sprintf('f_%d',i));
end
xlabel('log10(l)'); ylabel('Number of f evaluations');
title('Fibonacci: f-evals vs l');
legend('Location','best'); grid on; hold off;

%% === 6. Dichotomy with derivative: show typical run and counts ===
figure('Name','Dichotomy using derivative: endpoint trace and counts');
for i=1:3
    % use analytic derivative if available (df{i}), else pass empty to use numeric
    [af,bf,it,fcount,history] = dichotomy_derivative(f{i}, df{i}, a0, b0, default_l, default_eps);
    subplot(3,1,i);
    k=0:size(history,1)-1;
    plot(k, history(:,1), '-o', k, history(:,2), '-x');
    title(sprintf('Dichotomy w/ derivative f_%d: iterations=%d, f-evals=%d', i, it, fcount));
    xlabel('iteration k'); ylabel('interval endpoints'); legend('a_k','b_k','Location','best'); grid on;
end

%% === 7. Save results if you want ===
% save('work1_results.mat');

%% === Subfunctions implementing methods ===

function [a,b,iter,fcount,history] = dichotomy(f, a, b, l, eps)
% DICHOTOMY Simple interval halving search (no derivative)
% Inputs: f - function handle
%         a,b - initial interval endpoints
%         l - desired final interval length
%         eps - small positive number (distance from midpoint)
% Outputs: a,b final interval; iter number of iterations; fcount num f evals;
%         history - matrix of [a b] after each iteration (row per iter)

    if nargin < 5, eps = 1e-3; end
    iter = 0;
    fcount = 0;
    history = [];
    while (b - a) > l
        m = (a + b) / 2;
        x1 = m - eps;
        x2 = m + eps;
        % evaluate f (count these)
        f1 = f(x1); f2 = f(x2);
        fcount = fcount + 2;
        if f1 < f2
            % minimum is in [a, x2]
            b = x2;
        else
            % minimum is in [x1, b]
            a = x1;
        end
        iter = iter + 1;
        history(iter, :) = [a, b]; %#ok<AGROW>
        % safety: break if too many iterations
        if iter > 1e6, warning('dichotomy: too many iter'); break; end
    end
end

function [a,b,iter,fcount,history] = golden(f, a, b, l)
% GOLDEN Golden-section search (derivative-free) minimizing f on [a,b]
    phi = (sqrt(5) - 1) / 2; % ~0.618...
    % initial interior points
    x1 = b - phi * (b - a);
    x2 = a + phi * (b - a);
    f1 = f(x1); f2 = f(x2);
    fcount = 2; iter = 0; history = [];
    while (b - a) > l
        if f1 < f2
            b = x2;
            % move x2 <- x1, compute new x1
            x2 = x1;
            f2 = f1;
            x1 = b - phi * (b - a);
            f1 = f(x1);
            fcount = fcount + 1;
        else
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi * (b - a);
            f2 = f(x2);
            fcount = fcount + 1;
        end
        iter = iter + 1;
        history(iter, :) = [a, b]; %#ok<AGROW>
        if iter > 1e6, warning('golden: too many iter'); break; end
    end
end

function [a,b,iter,fcount,history] = fibonacci(f, a, b, l)
% FIBONACCI Fibonacci search for interval reduction to length <= l
    % compute minimal N such that F_N >= (b-a)/l
    fibs = [1,1]; % F1=1, F2=1
    while fibs(end) < (b - a) / l
        fibs(end+1) = fibs(end) + fibs(end-1); %#ok<AGROW>
        if length(fibs) > 1e4, error('fibonacci: too many fibs'); end
    end
    N = length(fibs);
    % We will use F_{N}, F_{N-1}, ...
    % initial two interior points
    k = 1;
    x1 = a + (fibs(N-2)/fibs(N)) * (b - a);
    x2 = a + (fibs(N-1)/fibs(N)) * (b - a);
    f1 = f(x1); f2 = f(x2);
    fcount = 2; iter = 0; history = [];
    % iterate N-2 times
    for k = 1:(N-2)
        if f1 < f2
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (fibs(N-k-2)/fibs(N-k)) * (b - a);
            f1 = f(x1);
            fcount = fcount + 1;
        else
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + (fibs(N-k-1)/fibs(N-k)) * (b - a);
            f2 = f(x2);
            fcount = fcount + 1;
        end
        iter = iter + 1;
        history(iter, :) = [a, b]; %#ok<AGROW>
    end
    % Final shrink if needed
    % Common implementations may perform an extra check to ensure interval <= l
end

function [a,b,iter,fcount,history] = dichotomy_derivative(f, df, a, b, l, eps)
% DICHOTOMY_DERIVATIVE Use derivative at midpoint to reduce interval
% If analytic df provided, we use it (no extra f evaluations).
% If df empty, use central difference (two f-evals) and count them.
    if nargin < 6, eps = 1e-3; end
    iter = 0; fcount = 0; history = [];
    while (b - a) > l
        m = (a + b) / 2;
        if isempty(df)
            % numerical derivative via central difference (counts two f-evals)
            h = 1e-6;
            fm1 = f(m - h); fm2 = f(m + h);
            fcount = fcount + 2;
            d = (fm2 - fm1) / (2*h);
        else
            % analytic derivative provided (no extra f evaluations counted)
            d = df(m);
        end
        if d > 0
            % derivative positive -> function increasing at m -> min to the left
            b = m;
        elseif d < 0
            a = m;
        else
            % stationary point found (rare numerically). shrink slightly.
            a = m - eps; b = m + eps;
        end
        iter = iter + 1;
        history(iter, :) = [a, b]; %#ok<AGROW>
        if iter > 1e6, warning('dichotomy_derivative: too many iter'); break; end
    end
end
