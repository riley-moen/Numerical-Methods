function u = CTCS(c, h, k, x0, xend, t0, tend, ic)

x = x0:h:xend;
t = t0:k:tend;
lambda = k^2/h^2;

A = zeros(length(x), length(x));
u = zeros(length(t), length(x));
A(1, 1) = 2;
for i = 2:length(x) - 1
    A(i, i - 1) = c^2*lambda;
    A(i, i) = 2 - 2*c^2*lambda;
    A(i, i + 1) = c^2*lambda;
end
A(end, end - 1) = 2;

u(1, :) = feval(ic, x);
u(2, :) = 0.5*A*u(1, :)';

for i = 3:length(t)
    u(i, :) = A*u(i - 1, :)' - u(i - 2, :)';
end

end