function u = CSCS(h, k, x0, xend, y0, yend, bc)

lambda = k^2 / h^2;
x = x0:h:xend;
y = y0:k:yend;
y1 = feval(bc, x);
x1 = feval(bc, y);

xsteps = length(x);
ysteps = length(y);
n = xsteps*ysteps;

b = zeros(xsteps, ysteps);
for i = 1:ysteps
    b(end, i) = x1(i);
end
for i = 1:xsteps
    b(i, end) = y1(i);
end
b = reshape(b, n, 1);

A = zeros(n, n);

for i = 1:xsteps
    A(i, i) = 1;
    A(end - i + 1, end - i + 1) = 1;
end

for i = (xsteps + 1):(n - xsteps)
    if b(i) ~= 0
        A(i, i) = 1;
    else
        A(i, i - 1) = -lambda;
        A(i, i) = 2 + 2*lambda;
        A(i, i + 1) = -lambda;
    end
end
A
u = A \ b;
u = reshape(u, xsteps, ysteps);
u = A
end
