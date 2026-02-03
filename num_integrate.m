function [T,M,S,G] = num_integrate(a, b, func, n)

% This function performs numerical integration on the values of func on the
% domain a to b, with n steps. 

h = (b-a)./n;

T = 0;
M = 0;
S = 0;
G = 0;
%Set up gaussian
[x, w] = gauss(n);
%Remap t to x
t = (b - a)/2.*x + (b+a)/2;

for i = 1:n
    %Adds the area under the trapezoidal rule
    T = T + feval(func, a + (i-1)*h) + feval(func, a+i*h);
    %Adds the area under each midpoint rectangle
    M = M + feval(func, a + (2*i-1)*h/2);
    if mod(i, 2) == 0
        %only takes the even points and applies simpsons quadrature to it
        S = S + feval(func, a + (i - 2)*h) + 4*feval(func, a + (i - 1)*h) + feval(func, a + i*h);
    end
    %Adds the guassian quadrature based on weights
    G = G + feval(func, t(i)).*w(i);
end

T = T * (h/2);
M = M * h;
S = S * (h/3);
G = G * (b-a)./2;

end