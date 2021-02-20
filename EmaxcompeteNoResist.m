function [dndt] = EmaxcompeteNoResist(t, y, param)

hill = param.hill; % Hill parameter
rMax = param.rMax; % Growth rate (hour^-1)
K = param.K; % Carrying capacity
rMin = param.rMin; % Maximum killing rate
epsiR = param.epsiR; % Growth rate depedence of kill rate 
gammaR = param.gammaR; % Growth rate for half kill rate % 0.1
w = param.w; % Antibiotic molecular weight
epsiK = param.epsiK; % Density-dependence of kill rate
kc = param.kc; % k0 in manuscript. Antibiotic concentration at which killing is half of max 
xi = param.xi; % Max increase density-dependent kill
gammaK = param.gammaK; % Density for k at half of its maximum
deltaB = param.deltaB; % Basal death rate (hour^-1)
deltaA = param.deltaA; % Antibiotic decay rate (hour^-1) % 0.05
lambda = param.lambda; %Phage production rate (hour^-1)
deltaV = param.deltaV; %Phage decay rate (hour^-1)
kD = param.kD; % Binding dissociation factor
exponent = param.exponent;

theta = param.theta; % Metabolic cost
%pDeath = param.pDeath; % Additional mortality from phage production 
phi = param.phi; % Antibiotic sequestration factor

Na = 6*10^23;

% y1 is phage-positive, y2 is phage-negative

B_p = y(1);
B_n = y(2);
V = y(3);
A = y(4);

N = B_p + B_n;

rmaxEff1 = rMax*(1 - N/K)^exponent*(1 - theta); 
rmaxEff2 = rMax*(1 - N/K)^exponent;

% Antibiotic concentration
Am = A*Na/w; 
Aeff = A - w*((Am + phi*V + kD) - ((Am + phi*V + kD)^2 - 4*Am*phi*V)^(1/2))/(2*Na);

% KC50
kEff1 = kc + epsiK*xi*(N/(N + gammaK));
kEff2 = kc + epsiK*xi*(N/(N + gammaK));

% Max kill rate
rMinEff1 = rMin*(1 - epsiR) + epsiR*rMin*(rmaxEff1/(rmaxEff1 + gammaR));
rMinEff2 = rMin*(1 - epsiR) + epsiR*rMin*(rmaxEff2/(rmaxEff2 + gammaR));

% Net growth rate
r1 = rmaxEff1 - (deltaB - rMinEff1/(1 + (kEff1/Aeff)^hill));
r2 = rmaxEff2 - (deltaB - rMinEff2/(1 + (kEff2/Aeff)^hill));

dB_p = r1*B_p;
if (B_p < 1) 
    dB_p = 0;
end

dB_n = r2*B_n;
if (B_n < 1) 
    dB_n = 0;
end

dV = dB_p*lambda - V*deltaV;

dA = - deltaA*A;

dndt = [dB_p;dB_n;dV;dA];

end