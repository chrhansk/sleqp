% Simple problem with single feasible point [0; 0]

function f = obj_val (x)
  f = x(1);
end

function g = obj_grad (x)
  g = [ 1
        0 ];
end

function c = cons_val (x)
  c = [ x(1) + x(2); x(1) - x(2) ];
end

function J = cons_jac (x)
  J = sparse([ 1 1; 1 -1 ]);
end

x_init          = [5; 5];  % The starting point.

funcs.obj_val        = @obj_val;
funcs.cons_val       = @cons_val;
funcs.obj_grad       = @obj_grad;
funcs.cons_jac       = @cons_jac;

options.cons_lb   = [0 0];   % Lower bounds on the constraint functions.
options.cons_ub   = [0 0];   % Upper bounds on the constraint functions.

% Enable damped BFGS
options.hess_eval   = sleqp.HessEval.DampedBFGS();

options.deriv_check = bitor(sleqp.DerivCheck.FirstObj(), sleqp.DerivCheck.FirstCons());

disp(options.deriv_check);

[x_act info] = sleqp.solve(x_init, funcs, options);

assert(info.status == sleqp.Status.Optimal());

assert(all(info.working_vars == sleqp.ActiveState.Inactive()));

assert(all(info.working_cons == bitor(sleqp.ActiveState.Lower(), sleqp.ActiveState.Upper())));
