clear; clc; close all;

K = 10; 

% Seed
rng(1888); 

u1_range = [-1, 2];
u2_range = [-2, 1];

% GA Parameters
pop_size = 100;
max_gen = 600;        
mutation_rate = 0.15;
selection_pressure = 4; 

% Training Data
num_train = 300;
u1_train = (u1_range(2)-u1_range(1)) * rand(num_train, 1) + u1_range(1);
u2_train = (u2_range(2)-u2_range(1)) * rand(num_train, 1) + u2_range(1);
y_train = target_func(u1_train, u2_train);

% Validation Data
num_val = 200;
u1_val = (u1_range(2)-u1_range(1)) * rand(num_val, 1) + u1_range(1);
u2_val = (u2_range(2)-u2_range(1)) * rand(num_val, 1) + u2_range(1);
y_val = target_func(u1_val, u2_val);

% Gene structure: [w, c1, sigma1, c2, sigma2]
num_genes = K * 5;

% Init
pop = zeros(pop_size, num_genes);
for i = 1:pop_size
    for k = 1:K
        idx = (k-1)*5;
        pop(i, idx+1) = randn(); % w
        pop(i, idx+2) = (u1_range(2)-u1_range(1))*rand() + u1_range(1); % c1
        pop(i, idx+3) = rand() + 0.2; % sigma1 (avoid too narrow initially)
        pop(i, idx+4) = (u2_range(2)-u2_range(1))*rand() + u2_range(1); % c2
        pop(i, idx+5) = rand() + 0.2; % sigma2
    end
end

best_fitness_history = [];
best_chrom = [];

fprintf('--- Starting GA optimization (Pop: %d, K: %d) ---\n', pop_size, K);

for gen = 1:max_gen
    % Evaluation (MSE)
    fitness = zeros(pop_size, 1);
    for i = 1:pop_size
        y_pred = evaluate_model(pop(i,:), u1_train, u2_train, K);
        fitness(i) = mean((y_train - y_pred).^2);
    end
    
    % Store best
    [min_fit, best_idx] = min(fitness);
    best_fitness_history = [best_fitness_history; min_fit];
    best_chrom = pop(best_idx, :);
    
    if mod(gen, 50) == 0
        fprintf('Generation %d | Best MSE: %.6f\n', gen, min_fit);
    end
    
    % Keep best
    new_pop = zeros(size(pop));
    new_pop(1, :) = best_chrom;
    
    % 2. Selection, Crossover, Mutation
    for i = 2:pop_size
        % Tournament Selection
        p1 = tournament_select(pop, fitness, selection_pressure);
        p2 = tournament_select(pop, fitness, selection_pressure);
        
        % Crossover
        alpha = rand(); 
        child = alpha * p1 + (1 - alpha) * p2;
        
        % Mutation
        if rand() < mutation_rate
            m_idx = randi(num_genes);
            % Add Gaussian noise
            child(m_idx) = child(m_idx) + 0.2 * randn(); 
        end
        
        % Constraints
        for k = 1:K
            base = (k-1)*5;
            
            child(base+2) = max(u1_range(1), min(u1_range(2), child(base+2)));
            
            % Clamp sigma1 (> 0.05 to prevent division by zero/spikes)
            child(base+3) = max(0.05, abs(child(base+3)));
            
            child(base+4) = max(u2_range(1), min(u2_range(2), child(base+4)));
            
            % Clamp sigma2
            child(base+5) = max(0.05, abs(child(base+5)));
        end
        
        new_pop(i, :) = child;
    end
    
    pop = new_pop;
end

% --- EVALUATION & PLOTTING ---
y_val_pred = evaluate_model(best_chrom, u1_val, u2_val, K);
val_mse = mean((y_val - y_val_pred).^2);

fprintf('\n');
fprintf('Final Training MSE: %.6f\n', min_fit);
fprintf('Validation MSE:     %.6f\n', val_mse);
fprintf('\n');


fprintf('\nPROPOSED ANALYTICAL EXPRESSION PARAMS (K=%d)\n', K);
fprintf('%-6s %-10s %-10s %-10s %-10s %-10s\n', 'Term', 'Weight', 'c1', 'sigma1', 'c2', 'sigma2');
for k = 1:K
    idx = (k-1)*5;
    w = best_chrom(idx+1);
    c1 = best_chrom(idx+2); s1 = best_chrom(idx+3);
    c2 = best_chrom(idx+4); s2 = best_chrom(idx+5);
    fprintf('%-6d %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n', k, w, c1, s1, c2, s2);
end
fprintf('--------------------------------------------------\n');

% Convergence
figure('Name', 'Convergence');
plot(best_fitness_history, 'LineWidth', 2);
xlabel('Generation'); ylabel('MSE');
title(['Convergence (Pop: ', num2str(pop_size), ', K: ', num2str(K), ')']);
grid on;
saveas(gcf, 'convergence_plot.png');


step = 0.05;
[U1_grid, U2_grid] = meshgrid(u1_range(1):step:u1_range(2), u2_range(1):step:u2_range(2));
Y_true_grid = target_func(U1_grid, U2_grid);

Y_pred_grid = zeros(size(U1_grid));
for k = 1:K
    idx = (k-1)*5;
    w = best_chrom(idx+1);
    c1 = best_chrom(idx+2); s1 = best_chrom(idx+3);
    c2 = best_chrom(idx+4); s2 = best_chrom(idx+5);
    exponent = - ( ((U1_grid - c1).^2)./(2*s1^2) + ((U2_grid - c2).^2)./(2*s2^2) );
    Y_pred_grid = Y_pred_grid + w .* exp(exponent);
end

% Comparison Surface
figure('Name', 'Comparison', 'Position', [100, 100, 1000, 400]);
subplot(1,2,1);
surf(U1_grid, U2_grid, Y_true_grid, 'EdgeColor', 'none');
title('True Function Surface');
xlabel('u1'); ylabel('u2'); zlabel('y');
view(45, 30); colorbar; shading interp;

subplot(1,2,2);
surf(U1_grid, U2_grid, Y_pred_grid, 'EdgeColor', 'none');
title(['Approximated (MSE: ' num2str(val_mse, '%.4f') ')']);
xlabel('u1'); ylabel('u2'); zlabel('y');
view(45, 30); colorbar; shading interp;
saveas(gcf, 'comparison_surface.png');

% Absolute Error
AbsError = abs(Y_true_grid - Y_pred_grid);
figure('Name', 'Error Distribution');
surf(U1_grid, U2_grid, AbsError, 'EdgeColor', 'none');
title(['Abs Error (|True - Pred|) | Max: ' num2str(max(max(AbsError)), '%.4f')]);
xlabel('u1'); ylabel('u2'); zlabel('Abs Error');
colormap(jet); colorbar; shading interp;
view(45, 30);
saveas(gcf, 'error_surface.png');

% Helpers

function y = target_func(u1, u2)
    y = sin(u1 + u2) .* sin(u2.^2);
end

function y_pred = evaluate_model(chrom, u1, u2, K)
    y_pred = zeros(size(u1));
    for k = 1:K
        idx = (k-1)*5;
        w = chrom(idx+1);
        c1 = chrom(idx+2);
        s1 = chrom(idx+3);
        c2 = chrom(idx+4);
        s2 = chrom(idx+5);

        exponent = - ( ((u1 - c1).^2)./(2*s1^2) + ((u2 - c2).^2)./(2*s2^2) );
        y_pred = y_pred + w .* exp(exponent);
    end
end

function parent = tournament_select(pop, fitness, k)
    pop_n = size(pop, 1);
    candidates_idx = randi(pop_n, k, 1);
    candidates_fit = fitness(candidates_idx);
    [~, best_loc] = min(candidates_fit);
    parent = pop(candidates_idx(best_loc), :);
end