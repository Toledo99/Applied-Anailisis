function [alpha] = bisect(fun,xk,pk)
% Input:
%   fun - function handle: R^n --> R
%   xk - current point: vector in R^n
%   pk - search direction: vector in R^n
% Output:
%   alpha - step size that minimizes fun(xk + alpha*pk)

a = 0;
b = 1;
c = (a+b)/2;

funlin = @(alpha) fun(xk+alpha*pk);

EPS = 1e-12;
MAXITER = 1e4;
k = 0;

while (b-a) >= EPS && k <= MAXITER
    y = (a+c)/2;
    fc = funlin(c);
    fy = funlin(y);
    if fy <= fc
        b = c;
        c = y;
        fc = fy;
    else
        z = (b+c)/2;
        fz = funlin(z);
        if fc <= fz
            a = y;
            b = z;
        else
            a = c;
            c = z;
            fc = fz;
        end
    end
    k = k+1;
end
alpha = c;
