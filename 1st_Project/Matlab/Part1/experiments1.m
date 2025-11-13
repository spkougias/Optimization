clear; close all; clc;

% functions from assignment
f1 = @(x) 5^x + (2 - cos(x))^2;
f2 = @(x) (x - 1)^2 + exp(1)*sin(x + 3);
f3 = @(x) exp(1) - (sin(x - 2) - 2)^2;

fs = {f1, f2, f3};
fnames = {'f1','f2','f3'};
a0 = -1; b0 = 3;
%-----------------------------------------------------------
%% A: l fixed, e varies
l_fixed = 0.01;
e_vals = [1e-6,1e-5,1e-4,1e-3,2e-3,4e-3];

figure('Name','Dichotomy: evaluations vs e (constatn l=0.01)');
for i=1:2
    evals_e = zeros(size(e_vals));
    for j=1:length(e_vals)
        results = dichotomy_no_deriv(fs{i}, a0, b0, l_fixed, e_vals(j), 1e5);
        evals_e(j) = results.evals;
    end
    plot(e_vals, evals_e, '-o','DisplayName',fnames{i});
    hold on;
end

xlabel('e');
ylabel('Number of f(x) evaluations');
title('Dichotomy: evaluations vs e (constant l = 0.01)');
legend('Location','best');
grid on;

%% A.b extreme values of e
l_fixed = 0.01;
e_vals = [1e-6,1e-5,1e-4,1e-3,2e-3,4e-3,6e-3,8e-3,1e-2,5e-2,1e-1];

figure('Name','Dichotomy: evaluations vs e (constatn l=0.01)');
for i=1:2
    evals_e = zeros(size(e_vals));
    for j=1:length(e_vals)
        results = dichotomy_no_deriv(fs{i}, a0, b0, l_fixed, e_vals(j), 1e5);
        evals_e(j) = results.evals;
    end
    plot(e_vals, evals_e, '-o','DisplayName',fnames{i});
    hold on;
end

set(gca,'XScale','log'); %Needed to be log so that the results can be visible
xlabel('e');
ylabel('Number of f(x) evaluations');
title('Dichotomy: evaluations vs e (constant l = 0.01)');
legend('Location','best');
grid on;
%-----------------------------------------------------------


%-----------------------------------------------------------
%% B: eps fixed, vary l
e_fixed = 1e-3;
l_vals = [1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1];

figure('Name','Dichotomy: evals vs l (eps=1e-3)');
for i=1:3
    evals_l = zeros(size(l_vals));
    for j=1:length(l_vals)
        results = dichotomy_no_deriv(fs{i}, a0, b0, l_vals(j), e_fixed, 1e5);
        evals_l(j) = results.evals;
    end
    plot(l_vals, evals_l, '-o','DisplayName',fnames{i}); hold on;
end
set(gca,'XScale','log'); set(gca,'YScale','log');
xlabel('l (final interval length)'); ylabel('Number of f(x) evaluations');
title('Dichotomy: evals vs l (\epsilon = 10^{-3})');
legend('Location','best'); grid on;
%-----------------------------------------------------------


%-----------------------------------------------------------
%% C: plot a_k and b_k vs k for several l values (for each function)
l_choices = [0.5, 0.1, 0.01];
for i=1:3
    figure('Name',['Interval endpoints vs k for ' fnames{i}]);
    for j=1:length(l_choices)
        results = dichotomy_no_deriv(fs{i}, a0, b0, l_choices(j), 1e-3, 1e5);
        k = 0:results.iter;
        subplot(length(l_choices),1,j);
        plot(k, results.a_list, '-o','DisplayName','a_k'); hold on;
        plot(k, results.b_list, '-s','DisplayName','b_k'); 
        xlabel('iteration k'); ylabel('endpoint value');
        title(sprintf('%s : l = %g   (iters = %d)', fnames{i}, l_choices(j), results.iter));
        legend('Location','best'); grid on;
    end
end
