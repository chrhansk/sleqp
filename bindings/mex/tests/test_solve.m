% HS71 test problem from the CUTest suite

function f = obj_val (x)
  f = x(1)*x(4)*sum(x(1:3)) + x(3);
end

function g = obj_grad (x)
  g = [ x(1)*x(4) + x(4)*sum(x(1:3))
        x(1)*x(4)
        x(1)*x(4) + 1
        x(1)*sum(x(1:3)) ];
end

function c = cons_val (x)
  c = [ prod(x); sum(x.^2) ];
end

function J = cons_jac (x)
  J = sparse([ prod(x)./x; 2*x ]); 
end

function H = hess (x, obj_dual, cons_dual)

  H = obj_dual*[ 2*x(4)             0      0   0;
                 x(4)               0      0   0;
                 x(4)               0      0   0;
                 2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ];

  H = H + cons_dual(1)*[    0          0         0         0;
                         x(3)*x(4)     0         0         0;
                         x(2)*x(4) x(1)*x(4)     0         0;
                         x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ];

  H = H + cons_dual(2)*diag([2 2 2 2]);
  H = sparse(H);
end

x0              = [1 5 5 1];  % The starting point.
options.var_lb  = [1 1 1 1];  % Lower bound on the variables.
options.var_ub  = [5 5 5 5];  % Upper bound on the variables.
options.cons_lb = [25  40];   % Lower bounds on the constraint functions.
options.cons_ub = [inf 40];   % Upper bounds on the constraint functions.

% Callbacks
funcs.objective         = @obj_val;
funcs.constraints       = @cons_val;
funcs.gradient          = @obj_grad;
funcs.jacobian          = @cons_jac;
funcs.hessian           = @hess;

% Run sleqp.
[x info] = sleqp.solve(x0, funcs, options);

disp("x:");
disp(x);
disp("info:");
disp(info);
