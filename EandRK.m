function [Euler, RK] = EandRK(f, h, t0, y0, finalt)

% This function solves the given differential equation using both Euler and
% Runge-Kutta 4 algorithms. 

steps = floor((finalt - t0) / h);
Euler = zeros(steps, 1);
Euler(1) = y0 + h*feval(f, t0, y0);
index = 1;
i = t0;

for j = 2:steps
    index = index + 1;
    Euler(index) = Euler(index - 1) + h*feval(f, i, Euler(index - 1));
    i = i + h;
end

RK = zeros(steps + 1, 1);
k = zeros(steps, 4);
index = 2;
RK(1) = y0;
k(index, 1) = feval(f, t0, y0);
k(index, 2) = feval(f, t0 + h/2, y0 + (h/2)*k(index, 1));
k(index, 3) = feval(f, t0 + h/2, y0 + (h/2)*k(index, 2));
k(index, 4) = feval(f, t0 + h, y0 + h*k(index, 3));

RK(index) = y0 + h*((1/6)*k(index, 1) + (1/3)*k(index, 2) + (1/3)*k(index, 3) + (1/6)*k(index, 4));
i = t0;

for j = 2:steps
    i = i + h;
    index = index + 1;
    k(index, 1) = feval(f, i, RK(index - 1));
    k(index, 2) = feval(f, i + h/2, RK(index - 1) + (h/2)*k(index, 1));
    k(index, 3) = feval(f, i + h/2, RK(index - 1) + (h/2)*k(index, 2));
    k(index, 4) = feval(f, i + h, RK(index - 1) + h*k(index, 3));
    RK(index) = RK(index - 1) + h*((1/6)*k(index, 1) + (1/3)*k(index, 2) + (1/3)*k(index, 3) + (1/6)*k(index, 4));
end

%RK = RK(end);
end


