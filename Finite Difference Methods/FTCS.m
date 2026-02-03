function u = FTCS(a, k, h, x0, xend, t0, tend, ic)

x = x0:h:xend;
t = t0:k:tend;
lambda = k / h;
v0 = feval(ic, x)';
A = zeros(length(x), length(x));
A(1, 1) = 1;
for i = 2:length(x) - 1
    A(i, i - 1) = a*lambda/2;
    A(i, i) = 1 - a*lambda/2;
    A(i, i + 1) = a * lambda/2;
end
A(end, end - 1) = 2*a*lambda;
A(end, end) = 1 - a*lambda;

u = zeros(length(t), length(x));
u(1, :) = v0;

for i = 2:length(t)
    v0 = A*v0;
    u(i, :) = v0';
end