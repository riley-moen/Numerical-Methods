clear all;

set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',20);

% This script uses an Ising model to simulate the magnetiziation of a
% metal. It compares many metrics of the model such as magnetization,
% specific heat, energy, magnetic susceptibility, and temperature.

Temps = linspace(0.5, 6.5, 50);

L = 15;
H = 0.1;
J = 1;

adj = adjMatrix(L);

EnCalc = @(s) -H*sum(s) - s'*(0.5*J*adj)*s;

EnArr = zeros([1 length(Temps)]);
MagArr = EnArr; E2Arr = EnArr; Mag2Arr = EnArr;

pass = 50;

for j = 1:length(Temps)

    T = Temps(j);

    lat = metropolis(ones([1 L^2]), 2*L, T, adj);
    
    Energy = 0;
    E2 = 0;
    Mag = 0;
    Mag2 = 0;

    for i = 1:pass
        U = metropolis(lat, L, T, adj);

        Esig = EnCalc(U');

        Energy = Energy + Esig;
        E2 = E2 + E2^2;

        Mag = Mag + (1/L^2)*sum(U);
        Mag2 = Mag + ((1/L^2)*sum(U))^2;

        lat = U;
    end
    EnArr(j) = Energy/pass;
    E2Arr(j) = (E2 - Energy^2)*(1/T^2);
    MagArr(j) = Mag/pass;
    Mag2Arr(j) = (Mag2 - Mag^2)*(1/T);
end

figure(1);
tiledlayout(1, 2)
nexttile
plot(Temps, EnArr)
title('Energy vs Temperature')
xlabel('Temperature')
ylabel('Energy')
nexttile
plot(Temps, E2Arr)
title('Specific Heat vs Temperature')
xlabel('Temperature')
ylabel('Specific Heat')

figure(2);
tiledlayout(1, 2)
nexttile
plot(Temps, MagArr)
title('Magnetization vs Temperature')
xlabel('Temperature')
ylabel('Magnetization')
nexttile
plot(Temps, Mag2Arr)
title('Suceptibility vs Temperature')
xlabel('Temperature')
ylabel('Susceptibility')


function A = adjMatrix(n)
    T = diag(ones(n-1,1), 1) + diag(ones(n-1,1), -1);
    T(1,n) = 1; T(n,1) = 1;
    
    I = speye(n);
    A = kron(T, I) + kron(I, T);
end

function sigma = metropolis(lat, L, T, adj)

    U = zeros([L length(lat)]);

    U(1, :) = lat;

    for i = 2:L^2
        k = random('unid', length(lat));

        lat0 = lat;
        
        lat(k) = lat(k)*-1;
    
        % de = 2*H*lat(k) + 2*lat(k)*(sum(lat) - lat(k));

        de = computeEnergy(lat, adj) - computeEnergy(lat0, adj);
    
        if de < 0
            U(i, :) = lat;
        elseif random('Uniform', 0, 1) <= exp(-de/T)
            U(i, :) = lat;
        else
            lat(k) = lat(k)*-1;
        end
    end
    
    U = U(any(U, 2), :);
    sigma = U(end, :);
end

function energy = computeEnergy(s, adj)
    H = 0.1;
    J = 1;
    energy = -H*sum(s) - s*(0.5*J*adj)*s';
end