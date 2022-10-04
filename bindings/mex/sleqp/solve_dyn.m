function [x,info] = solve_dyn(x0, funcs, options)
  [x info] = sleqp.extension('solve_dyn', x0, funcs, options);
end
