function [intervals, eval_num] = golden_section_method(f, a_initial, b_initial, l)
 
    % --- Setup ---
    intervals = [a_initial, b_initial]; % history of [a, b]
    gamma = 0.618034;
    x1 = a_initial + (1 - gamma)*(b_initial - a_initial);
    x2 = a_initial + gamma*(b_initial - a_initial);
    f_x1 = f(x1);
    f_x2 = f(x2);
    eval_num = 2; % The first two evaluations
    
    % Loop until less than l
    while ~(intervals(end,2) - intervals(end,1) < l)
        if f_x1 > f_x2
            % New interval [x1, b]
            intervals = [intervals; x1, intervals(end,2)];
            
            x1 = x2;
            f_x1 = f_x2;
            
            x2 = intervals(end,1) + gamma*(intervals(end,2) - intervals(end,1));
            f_x2 = f(x2);
            
        elseif f_x1 < f_x2
            % New interval [a, x2]
            intervals = [intervals; intervals(end,1), x2];

            x2 = x1;
            f_x2 = f_x1;

            x1 = intervals(end,1) + (1 - gamma)*(intervals(end,2) - intervals(end,1));
            f_x1 = f(x1);
        end
        eval_num = eval_num + 1; % One evaluation per loop
    end
end