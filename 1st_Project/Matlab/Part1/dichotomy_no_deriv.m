% Inputs:
%   f       - function fi(x)
%   a,b     - [a,b]
%   l       - b-a < l
%   eps     - (a+b)/2 +- e
%   max - maximum number of iterations
%
% Outputs:
%   .a_hist, .b_hist  - vectors of a and b potition at each iteration
%   .xmid             - final midpoint
%   .fmin             - f at final midpoint
%   .iter, .evals     - number of iterations and total evaluations
function result = dichotomy_no_deriv(f,a,b,l,e,maxi)
    
    
    if nargin<6, maxi = 1e5; end         % Max iterations set to 10000 if not given
    
    a_list = a; b_list = b;
    iter = 0;     % Iterations
    evals = 0;    % Number of evaluations
    
    if e <= 0
        error('eps must be > 0');
    end
    
    while (b - a) > l && iter < maxi
        mid = (a + b)/2;
        x1 = mid - e;
        x2 = mid + e;
    
        % make sure that x1,x2 are inside [a,b]
        x1 = max(x1,a); x2 = min(x2,b);
    
        f1 = f(x1); f2 = f(x2);
        evals = evals + 2;
    
        if f1 <= f2
            b = x2;
        else
            a = x1;
        end
        iter = iter + 1;
        % vectors of a and b positions in every iteration
        a_list(end+1) = a; 
        b_list(end+1) = b;
    end
    
    xmid = (a + b)/2;
    fmin = f(xmid); evals = evals + 1;
    
    result.a_list = a_list;
    result.b_list = b_list;
    result.xmid = xmid;
    result.fmin = fmin;
    result.iter = iter;
    result.evals = evals;
    result.final_interval = b - a;
end
