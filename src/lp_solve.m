function [x, ok] = lp_solve(c, A, b, Aeq, beq)
%LP_SOLVE Portable LP wrapper:  min c'x  s.t.  A x <= b,  Aeq x = beq.
%   Uses linprog (MATLAB Optimization Toolbox) if available, otherwise
%   GLPK (built into GNU Octave). Variables are free (unbounded).
n = numel(c);
if exist('linprog', 'file') == 2 && exist('optimoptions', 'file') == 2
    o = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');
    [x, ~, flag] = linprog(c, A, b, Aeq, beq, [], [], o);
    ok = (flag == 1);
elseif exist('glpk', 'file')
    AA    = [A; Aeq];
    bb    = [b; beq];
    ctype = [repmat('U', size(A,1), 1); repmat('S', size(Aeq,1), 1)].';
    vtype = repmat('C', 1, n);
    lb    = -Inf(n, 1);  ub = Inf(n, 1);
    [x, ~, err] = glpk(c, AA, bb, lb, ub, ctype, vtype, 1);
    ok = (err == 0);
else
    error('No LP solver found (need linprog or glpk).');
end
end
