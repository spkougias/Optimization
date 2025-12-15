function [x_history, k] = steepest_descent(f, vars, x_start, epsilon, gamma_initial, gamma_mode)
    % Implements Steepest Descent Optimization
    % Inputs:
    %   f: Symbolic function
    %   vars: Symbolic variables [x, y]
    %   x_start: Starting point [x0, y0]
    %   epsilon: Convergence tolerance
    %   gamma_initial: Initial step size value
    %   gamma_mode: 1=Fixed, 2=Armijo, 3=Exact(Bisection)
    % Output:
    %   x_history: Matrix of visited points
    %   k: Number of iterations performed

    syms gam_sym;
    
    % Gradient
    grad_f = gradient(f, vars);
    
    % Initialize
    x_current = x_start;
    x_history = x_start;
    k = 0;
    MAX_ITERS
    
    grad_val = double(subs(grad_f, vars, x_current));
    
    % Armijo Parameters
    alpha = 0.001; 
    beta = 0.5; 

    while norm(grad_val) > epsilon && k < MAX_ITERS
        
        d_k = -grad_val';

        % Step Size
        switch gamma_mode
            case 1 % Fixed Step
                gamma_k = gamma_initial;
                
            case 2 % Armijo Rule
                mk = 0;
                gamma_try = gamma_initial;
                while true
                    x_next_try = x_current + gamma_try * d_k;
                    
                    f_curr = double(subs(f, vars, x_current));
                    f_next = double(subs(f, vars, x_next_try));
                    
                    limit = f_curr + alpha * gamma_try * (d_k * grad_val);
                    
                    if f_next <= limit
                        gamma_k = gamma_try;
                        break;
                    end
                    mk = mk + 1;
                    gamma_try = gamma_initial * (beta^mk);
                    
                    if mk > 20, gamma_k = gamma_try; break; end % Safety
                end
            otherwise % Exact Search (Minimizing f along direction)
                f_gamma_sym = subs(f, vars, x_current + gam_sym * d_k);
                df_gamma_sym = diff(f_gamma_sym, gam_sym);
                df_gamma_handle = matlabFunction(df_gamma_sym);
                
                [intervals, ~] = bisection_with_dervs(df_gamma_handle, 0, gamma_initial, 0.001);
                gamma_k = (intervals(end, 1) + intervals(end, 2)) / 2;
        end
        
        % Update Position
        x_current = x_current + gamma_k * d_k;
        
        x_history = [x_history; x_current];
        grad_val = double(subs(grad_f, vars, x_current));
        k = k + 1;
    end
end