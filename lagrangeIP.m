function y = lagrangeIP(x, f, xi)

% This function uses Lagrange interpolation on the values of f to find a
% function y with xi number of values

n = length(x);
m = length(xi);
y = zeros(1, n);

for i = 1:n
    P = ones(m, 1);
    for j = 1:n
        if i ~= j
            P = P.*(xi' - x(j))/(x(i) - x(j));
        end
    end
    y = y + f(i).*P;
end
end