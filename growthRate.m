%%% This function calculates the bacterial growth rate.
% Arguments: Bacterial density 'B', Viral density 'V', Antibiotic concentration 'A' and parameter
% list 'param' and identity of strain 'Pf' (0 for Pf- and 1 for Pf+)
% Return: Growth rate 'r'

function [r] = growthRate(B,V,A,param,Pf)

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
pDeath = param.pDeath; % Additional mortality from phage production 
phi = param.phi; % Antibiotic sequestration factor
Na = 6*10^23;

% Effective growth rate
rmaxEff = rMax*(1-B/K)^exponent*(1-theta*Pf); 

% Effective antibiotic concentration
Am = A*Na/w; 
Aeff = A - w*((Am + phi*V + kD) - ((Am + phi*V + kD)^2 - 4*Am*phi*V)^(1/2))/(2*Na);

% Effective KC50 
kEff = kc + epsiK*xi*(B/(B + gammaK));

% Effective max killing rate
rMinEff = rMin*(1 - epsiR) + epsiR*rMin*(rmaxEff/(rmaxEff + gammaR));

% Net growth rate
r = rmaxEff - (deltaB - rMinEff/(1 + (kEff/Aeff)^hill));

