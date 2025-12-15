function [x_history, k] = steepest_descent(f, vars, x_start, epsilon, gamma, max_iters)
    % Implements Steepest Descent Optimization
    % Inputs:
    %   f: Symbolic function
    %   vars: Symbolic variables [x, y]
    %   x_start: Starting point [x0, y0]
    %   epsilon: Convergence tolerance
    %   gamma: Initial step size value
    %   max_iters: The maximum number of iterations
    % Output:
    %   x_history: Matrix of visited points
    %   k: Number of iterations performed

    grad_f = gradient(f, vars);
    
    x_current = x_start;
    x_history = x_start;
    k = 0;

    grad_val = double(subs(grad_f, vars, x_current));
    
    while norm(grad_val) > epsilon && k < max_iters

        d_k = -grad_val'; 
        
        % Update step
        x_next = x_current + gamma * d_k;
        
        % Update position and history
        x_history = [x_history; x_next];
        x_current = x_next;

        grad_val = double(subs(grad_f, vars, x_current));
        k = k + 1;
    end
    
    % Warning if not converge
    if k == max_iters
        fprintf('Error: max iterations (%d) reached for gamma = %.2f.\n', max_iters, gamma);
    end
end