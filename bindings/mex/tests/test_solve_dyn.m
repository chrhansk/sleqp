% Dynamic version of the Rosenbrock problem

function test_solve_dyn()
  x_init          = [0; 0];        % The starting point.
  options.var_lb  = [-inf; -inf];  % Lower bound on the variables.
  options.var_ub  = [ inf;  inf];  % Upper bound on the variables.

  % Callbacks
  funcs.eval             = @func_eval;
  funcs.obj_grad         = @obj_grad;
  funcs.hess             = @hess;

  % Set rng seed
  rand("seed", 42);

  global a b;

  a = 1;
  b = 100;

  % Run sleqp.
  [x_act info] = sleqp.solve_dyn(x_init, funcs, options);

  assert(info.status == sleqp.Status.Optimal());

  x_exp = [1.; 1.];

  diff = norm(x_exp - x_act, Inf);

  tolerance = 1e-6;

  assert(diff < tolerance);
end

function [obj, cons, error] = func_eval (x, error_info)
  global a b;

  noise = 2*rand() - 1.;
  factor = error_info.error_bound / error_info.obj_weight;
  error = factor * noise;

  orig_obj = (a - x(1))^2 + b*(x(2) - x(1)^2)^2;

  obj = orig_obj + error;

  error = abs(error);

  cons = [];
end

function g = obj_grad (x, error_info)
  global a b;

  g = [(4*b*x(1)*(x(1)^2 - x(2)) + 2*x(1) - 2*a)
       (-2*b*(x(1)^2 - x(2)))];
end

function H = hess (x, cons_dual, error_info)
  global a b;

  H = [(8*b*x(1)^2 + 4*b*(x(1)^2 - x(2)) + 2) 0;
       (-4*b*x(1)) (2*b)];

  H = sparse(H);
end
