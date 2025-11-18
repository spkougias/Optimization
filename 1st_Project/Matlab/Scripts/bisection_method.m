function [intervals, eval_num] = bisection_method(f, a_initial, b_initial, epsilon, l)

    intervals = [a_initial, b_initial]; % history of [a, b]
    eval_num = 0; % Number of Evaluations
    
    
    while ~(intervals(end,2) - intervals(end,1) < l)
        
        x1 = (intervals(end,1) + intervals(end,2))/2 - epsilon;
        x2 = (intervals(end,1) + intervals(end,2))/2 + epsilon;
        
        f_x1 = f(x1);
        f_x2 = f(x2);
        eval_num = eval_num + 2; % +2 Evaluations
        
        if f_x1 < f_x2
            intervals = [intervals; intervals(end,1), x2]; % Keep [a, x2]
        else
            intervals = [intervals; x1, intervals(end,2)]; % Keep [x1, b]
        end
    end
end