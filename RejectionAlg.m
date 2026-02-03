clear all;

set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',20);

% This script uses the Rejection algorithm to generate
% pseudo-random numbers from a given distribution. It then uses a
% Chi-Square test to see how close the generated numbers are to the given
% distribution.

f = @(x) (x.^4.* exp(-x))/24;

ginv = @(x) -5*log(x);

g = @(x) exp(-(0.2.*x));

p = @(x) (x.^4.* exp(-x))/24;


n = 100000;
x = [];

for i = 1:n
    y = ginv(random('Uniform', 0, 1));
    z = random('Uniform', 0, 1);

    if z <= (f(y)/g(y))
        x = [x; y];
    end
end


figure
set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',20)
histogram(x, 20,'Normalization', 'pdf');
hold on;
x = linspace(0, 20, 10000);
plot(x, f(x))
plot(x, g(x))
title('Rejection Method for $f(x) = \frac{x^4e^{-x}}{\Gamma (5)}')

figure
h = histogram(x, 20);

chisqr = 0;

for i = 1:h.NumBins
    r = linspace(h.BinEdges(i), h.BinEdges(i+1), 100);
    t = trapz(p(r));

    chisqr = chisqr + ((h.Values(i) - n*t)/(n*t)).^2;
end

alpha = 0.01;
a = 0.0001;

sig = sqrt((1/(n-1)*(sumsqr(x) - n*(mean(x))^2)));

minN = 2*((erf(alpha)*sig)/a)^2;