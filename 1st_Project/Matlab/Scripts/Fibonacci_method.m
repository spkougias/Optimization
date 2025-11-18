function [intervals, eval_num] = Fibonacci_method(f, a_initial, b_initial, epsilon, l)
    % l = specified accuracy
    intervals = [a_initial, b_initial]; % history of [a, b]
    
    % Fibonacci sequence
    fib = [1, 1];
    while fib(end) <= (b_initial - a_initial) / l
        fib = [fib, fib(end) + fib(end - 1)];
    end
    n = length(fib) - 1; % get index

    % Initial Points
    x1 = a_initial + (fib(n-2+1)/fib(n+1))*(b_initial - a_initial);
    x2 = a_initial + (fib(n-1+1)/fib(n+1))*(b_initial - a_initial);
    f_x1 = f(x1);
    f_x2 = f(x2);
    eval_num = 2; % The first two evaluations
    
    for k = 1:n-2
        if f_x1 > f_x2
            % New interval [x1, b]
            intervals = [intervals; x1, intervals(end,2)];
            x1 = x2;
            x2 = intervals(end,1) + (fib(n-k-1+1)/fib(n-k+1))*(intervals(end,2) - intervals(end,1));
            
            if k ~= n-2
                f_x1 = f_x2;
                f_x2 = f(x2);
            else
                x2 = x1 + epsilon;
                f_x1 = f_x2;
                f_x2 = f(x2);
            end
        elseif f_x1 < f_x2
            % New interval [a, x2]
            intervals = [intervals; intervals(end,1), x2];
            x2 = x1;
            x1 = intervals(end,1) + (fib(n-k-2+1)/fib(n-k+1))*(intervals(end,2) - intervals(end,1));
            
            if k ~= n-2
                f_x2 = f_x1;
                f_x1 = f(x1);
            else
                x2 = x1 + epsilon;
                f_x1 = f(x1);
                f_x2 = f(x2);
            end
        end
        eval_num = eval_num + 1; % One evaluation per loop
    end

    if f_x1 > f_x2
        intervals = [intervals; x1, intervals(end,2)];
    else
        intervals = [intervals; intervals(end,1), x2];
    end
end