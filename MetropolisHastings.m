clear all;

set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',20);


% This script uses the Metropolis-Hastings algorithm to generate
% pseudo-random numbers from a given distribution. It then uses a
% Chi-Square test to see how close the generated numbers are to the given
% distribution.

n = 100000;

x0 = 5;
x = [];
acc = [];
x(1) = x0;

f = @(x) 0.5*exp(-abs(x));

p = @(x) 0.5*exp(-abs(x));

delta = 4.5;

for i = 1:n

    mu = random('Uniform', -delta, delta);
    x1 = x0 + mu;

    alp = f(x1)/f(x0);

    if alp >= 1
        x = [x; x1];
        x0 = x1;
        acc = [acc; 1];
    elseif alp >= random('Uniform', 0, 1)
        x = [x; x1];
        x0 = x1;
        acc = [acc; 1];
    else
        x = [x; x0];
        acc = [acc; 0];
    end
end

[cor, lags] = autocorr(x, NumLags=100);

figure
tiledlayout(1, 2)
nexttile
histogram(x, 41,'Normalization', 'pdf');
hold on;
x = linspace(-5, 5, 10000);
plot(x, f(x))
xlim([-5 5])
title('Metropolis-Hastings Method for $p(x) = \frac{1}{2} \exp^{-|x|}$')
xlabel('x')
ylabel('Density')
legend('Generated Distribution', 'Target Distribution')
grid on;
nexttile
scatter(lags, cor, 'filled')
yline(0, 'LineWidth', 2, 'Color', 'r')
xlim([lags(1) lags(end)])
ylim([-0.2, 1])
title('Autocorrelation of Markov Chain')
ylabel('Autocorrelation')
xlabel('Lag')
grid on;

figure
per = zeros(n, 1);
for i = 1:n
    per(i) = sum(acc(1:i))/i;
end
plot(1:n, per)



figure
h = histogram(x, 40);

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