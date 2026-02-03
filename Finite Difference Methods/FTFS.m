function u = FTFS(a, k, h, x0, xend, t0, tend, ic)

x = x0:h:xend;
t = t0:k:tend;
lambda = k/h;

A = zeros(length(x), length(x));
A(1, 1) = 1;
for i = 2:length(x) - 1
    A(i, i) = 1 + a*lambda;
    A(i, i + 1) = -a*lambda;
end
A(end, end) = 1 + a*lambda;
u = zeros(length(t), length(x));
u(1, :) = feval(ic, x);
for i = 2:length(t)
    u(i, :) = A*u(i - 1, :)';
end