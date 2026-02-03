function u = FTBS(a, k, h, x0, xend, t0, tend, ic)

x = x0:h:xend;
t = t0:k:tend;

xlength = length(x);
tlength = length(t);

lambda = k/h;
A = zeros(xlength, xlength);
A(1, 1) = 1;

for i = 2:xlength
    A(i, i - 1) = a*lambda;
    A(i, i) = 1 - a*lambda;
end

u = zeros(tlength, xlength);
v0 = feval(ic, x);
u(1, :) = v0;

for i = 2:tlength
    v0 = A*v0';
    u(i, :) = v0;
end