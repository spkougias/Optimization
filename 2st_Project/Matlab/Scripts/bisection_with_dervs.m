function [intervals, eval_num] = bisection_with_dervs(df, a_initial, b_initial, l)

    intervals = [a_initial, b_initial]; % history of [a, b]
    eval_num = 0; % Number of Evaluations
    
    while ~(intervals(end,2) - intervals(end,1) < l)
        
        xk = (intervals(end,1) + intervals(end,2))/2; % midpoint
        
        df_xk = df(xk);
        eval_num = eval_num + 1; % +1 evaluation
        
        if df_xk == 0
            intervals = [intervals; xk, xk]; % Found minimum
            break;
        elseif df_xk > 0
            intervals = [intervals; intervals(end,1), xk]; % Keep left half
        else
            intervals = [intervals; xk, intervals(end,2)]; % Keep right half
        end
    end
end