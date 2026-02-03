function u = BTCS(a, h, k, x0, xend, t0, tend, ic)
% Computes the backwards time, center space finite difference method for a
% 1D space differential equation.

x = x0:h:xend;
t = t0:k:tend;
mu = k/h^2;

A = zeros(length(x), length(x));
A(1, 1) = 1 + 2*a^2*mu;
A(1, 2) = -a^2*mu;
A(end, end) = 1 + 2*a^2*mu;
A(end, end - 1) = -a^2*mu;
for i = 2:length(x) - 1
    A(i, i - 1) = -a^2*mu;
    A(i, i + 1) = -a^2*mu;
    A(i, i) = 1 + 2*a^2*mu;
end

u = zeros(length(t), length(x));
v0 = feval(ic, x)';
u(1, :) = v0;
for i = 2:length(t)
    v = A \ v0;
    u(i, :) = v;
    v0 = v;
end

