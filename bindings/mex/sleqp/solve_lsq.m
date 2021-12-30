function [x,info] = solve_lsq(x0, funcs, options)
  [x info] = sleqp._extension("solve_lsq", x0, funcs, options);
end
