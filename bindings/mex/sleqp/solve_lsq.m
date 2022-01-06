function [x,info] = solve_lsq(x0, funcs, options)
  [x info] = sleqp.extension('solve_lsq', x0, funcs, options);
end
