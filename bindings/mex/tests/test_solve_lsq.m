% LSQ variant of the Rosenbrock function


function r = lsq_res(x)
    a = 1;
    b = 100;
    r = [a-x(1); sqrt(b)*(x(2)-x(1)^2)];
end

function p = for_prod(x, d)
    a = 1;
    b = 100;     
    p = [-d(1); sqrt(b)*(-2*x(1)*d(1) + d(2))];
end

function p = adj_prod(x, d)
    a = 1;
    b = 100;         
    p = [-d(1) - 2*x(1)*sqrt(b)*d(2) sqrt(b)*d(2)];
end

x_init = [0; 0];  % The starting point.

options = struct;

% Callbacks
funcs.lsq_residuals   = @lsq_res;
funcs.lsq_jac_forward = @for_prod;
funcs.lsq_jac_adjoint = @adj_prod;

% Run sleqp.
[x_act info] = sleqp.solve_lsq(x_init, funcs, options);

assert(info.status == sleqp.Status.Optimal());

x_exp = [1; 1];

diff = norm(x_exp - x_act, Inf);

tolerance = 1e-8;

assert(diff < tolerance);

