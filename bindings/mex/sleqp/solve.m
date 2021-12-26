function [x,info] = solve(x0, funcs, options)
  [x info] = sleqp._extension("solve", x0, funcs, options);
end
