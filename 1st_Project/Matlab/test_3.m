

f{1} = @(x) 5 + (2 - cos(x)).^2;             % f1(x)
f{2} = @(x) (x - 1).^2 + exp(1).*sin(x + 3); % f2(x)
f{3} = @(x) exp( - (sin(x - 2) - 2) );       % f3(x)

% Inputs:

a0 = -1;            % - starting position
b0 = 3;             % - ending position
l = 0.01;           % - final required interval length (stop when b-a <= l)
eps  = 0.0001;      % - small separation around midpoint (epsilon > 0)
maxiter = 1000;     % - maximum number of iterations (safety)

res = dichotomy_no_deriv(f{1},a0,b0,l,0.0001,1000);
