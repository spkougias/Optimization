function [x_history, k] = projected_steepest_descent(f, vars, x_initial, epsilon, gamma, s, a, b, diminish_s)
    % Implements the projected steepest descent method with step sizes
    % gamma and s to find the minimum of a function f of
    % multiple variables vars over the region X = {x ∈ R^n : a_i ≤ x_i ≤ b_i, where i = 1, 2, ..., n},
    % and where a and b are the lower and upper bounds for each dimension.
    % Starts from an intitial point x_initial ∈ X and outputs all points
    % visited by the method until the norm of the gradient is smaller than
    % the specified tolerance epsilon.

    grad_f = gradient(f, vars);
    x_current = x_initial;

    % Project x_initial onto [a, b]
    x_current = min(max(x_current, a), b);
    
    x_history = x_current;
    
    % Thema 4, record the adjusted start
    if ~isequal(x_current, x_initial)
        x_history = [x_initial; x_current];
    end
    
    grad_val = double(subs(grad_f, vars, x_current));
    norm_grad = norm(grad_val);
    
    k = 0;
    max_iters = 500;
    
    while norm_grad >= epsilon && k < max_iters
        k = k + 1;
        
        z_k = x_current - s * grad_val';
        x_bar = min(max(z_k, a), b);    % project onto [a,b]
        
        % Update
        x_next = x_current + gamma * (x_bar - x_current);
        % Save history
        x_history = [x_history; x_next];
        x_current = x_next;

        grad_val = double(subs(grad_f, vars, x_current));
        norm_grad = norm(grad_val);
        
        % Thema 3, Diminishing s
        if diminish_s
            s = s * 0.9;
        end
    end
end