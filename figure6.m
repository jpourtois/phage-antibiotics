% Effect of initial conditions: Figure 6

%% Parameter values

% Define values for fixed parameters

param.hill = 0.8; % Hill parameter
param.rMax = 1; % Growth rate (hour^-1)
param.K = 10^8; % Carrying capacity
param.rMin = -12; % Maximum killing rate
param.epsiR = 0.5; % Growth rate depedence of kill rate 
param.gammaR = 0.2; % Growth rate for half kill rate % 0.1
param.w = 468*10^6; % Antibiotic molecular weight
param.epsiK = 1; % Density-dependence of kill rate
param.kc = 3.68; % k0 in manuscript. Antibiotic concentration at which killing is half of max 
param.xi = 10; % Max increase density-dependent kill
param.gammaK = 5*10^6; % Density for k at half of its maximum
param.deltaB = 0.05; % Basal death rate (hour^-1)
param.deltaA = 0.1; % Antibiotic decay rate (hour^-1) % 0.05
param.lambda = 1; %Phage production rate (hour^-1)
param.deltaV = 0.1; %Phage decay rate (hour^-1)
param.kD = 10^(14);% Binding dissociation factor 10^13-10^15
param.exponent = 1;

% Stochastic simulation? (0 for No, 1 for yes)

param.randomR = 0; % Random r?

% Default value for varying parameters 

param.theta = 0.2; % Metabolic cost
param.pDeath = 1; % Additional mortality from phage production 
param.phi = 10^6; % Antibiotic sequestration factor
param.aMax = 1; % Peak antibiotic concentration

param.MAX_T = 120; % Length of simulation(hour^-1)
MAX_T = param.MAX_T;

% Calculate MIC

mic = param.kc/(-param.rMin/(param.rMax - param.deltaB) - 1)^(1/param.hill);

%% Determine range for sensitivity analysis and generate random values

pMin(1) = 0.6; % Hill 
pMax(1) = 2.5;

pMin(2) = 0.5; % rMax
pMax(2) = 1.5;

pMin(3) = 2; % rMin
pMax(3) = 15;

pMin(4)= 0.1; % epsiR
pMax(4) = 0.9;

pMin(5)= 0.05; % gammaR
pMax(5) = 0.5;

pMin(6) = 10^(-1); % epsiK
pMax(6) = 1;

pMin(7) = 2; % xi
pMax(7) = 20;

pMin(8) = 10^6; % gammaK
pMax(8) = 10^7;

pMin(9) = 0.025; % deltaB
pMax(9) = 0.25;

pMin(10) = 0.05; % deltaA
pMax(10) = 0.5;

pMin(11) = 0.1; %lambda
pMax(11) = 10;

pMin(12) = 0.01; %deltaV
pMax(12) = 1;

pMin(13) = 10^(13); % kD
pMax(13) = 10^(15);

pMin(14) = 0.05;%theta
pMax(14) = 0.5;

pMin(15) = 10^5;%phi
pMax(15) = 10^7;

n_par = size(pMin,2);

%Log transform min and max
logMin = log(pMin);
logMax = log(pMax);

% Find range spread in log space
diffLog = logMax - logMin;

% Number of data sets generated
reps = 100;

% Specify seed for reproducibility
rng(2020);

% Generate 'reps'random sets of numbers, ranging from 0 to 1
randomUnscaled = lhsdesign(reps,n_par);

% Scale the random number with the range spread in log space
randomdiffLog = zeros(reps,n_par);
 
for n = 1:n_par
    randomdiffLog(:,n) = randomUnscaled(:,n)*diffLog(n);    
end

randomLogScaled = logMin + randomdiffLog;

%Get out of log space
randomScaled = exp(randomLogScaled);

%% Same initial conditions for Pf+ and Pf-, phages already present

range_dose = 10.^(-2:0.05:2);

[percent_win_1, percent_win_n_1, range_init_1, p_extinct_1] = fc_var_par_median(randomScaled, 'init', 1);

%%
[percent_win_0, range_init_0, p_extinct_0, p_extinct_p_0, p_extinct_n_0] = fc_var_par_median(randomScaled, 'init', 0);

%% Antibiotic sequestration by phages

% Define values for fixed parameters

param.hill = 0.8; % Hill parameter
param.rMax = 1; % Growth rate (hour^-1)
param.K = 10^8; % Carrying capacity
param.rMin = -12; % Maximum killing rate
param.epsiR = 0.5; % Growth rate depedence of kill rate 
param.gammaR = 0.2; % Growth rate for half kill rate % 0.1
param.w = 468*10^6; % Antibiotic molecular weight
param.epsiK = 1; % Density-dependence of kill rate
param.kc = 3.68; % k0 in manuscript. Antibiotic concentration at which killing is half of max 
param.xi = 10; % Max increase density-dependent kill
param.gammaK = 5*10^6; % Density for k at half of its maximum
param.deltaB = 0.05; % Basal death rate (hour^-1)
param.deltaA = 0.1; % Antibiotic decay rate (hour^-1) % 0.05
param.lambda = 1; %Phage production rate (hour^-1)
param.deltaV = 0.1; %Phage decay rate (hour^-1)
param.kD = 10^(14);% Binding dissociation factor 10^13-10^15
param.exponent = 1;

param.theta = 0.2; % Metabolic cost
param.pDeath = 1; % Additional mortality from phage production 
param.phi = 10^7; % Antibiotic sequestration factor
param.aMax = 1; % Peak antibiotic concentration

Na = 6*10^23;

v_range = 10.^(2:0.1:10);
A_range = 10.^(-2:0.05:2);

Aeff = zeros(size(v_range,2), size(A_range,2));


for i = 1:size(v_range,2)
    
    V = v_range(i);
    
    for j = 1:size(A_range,2)
        
        A = A_range(j);
        
        Am = A*Na/param.w; 
        Aeff(i,j) = A - param.w*((Am + param.phi*V + param.kD) - ((Am + param.phi*V + param.kD)^2 - 4*Am*param.phi*V)^(1/2))/(2*Na);
        Aeff(i,j) = 1 - Aeff(i,j)./A;
        
        
    end 
end

figure()
s = surface(A_range,v_range, Aeff, 'FaceAlpha',1)
colorbar
hold on
%shading interp
s.EdgeColor = 'none';
ylim([10^6 10^10])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%xlim([0.1 0.9])
xlabel('AMC (ug/ml)')
ylabel('Pf concentration (PFU/ml)')
title('Antibiotic sequestered (%)')



%%

figure;
subplot(1,2,1)
s = surface(range_init_1, range_dose./mic, percent_win_1 - percent_win_n_1 , 'FaceAlpha',1)
colorbar
hold on
shading interp
s.EdgeColor = 'none';
ylim([10^0 10^2])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%xlim([0.1 0.9])
xlabel('Initial density (CFU/ml)')
ylabel('AMC (xMIC)')
%title('Difference in runs (%)')

subplot(1,2,2)
s = surface(A_range,v_range, Aeff, 'FaceAlpha',1)
colorbar
hold on
%shading interp
s.EdgeColor = 'none';
ylim([10^6 10^10])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%xlim([0.1 0.9])
xlabel('AMC (ug/ml)')
ylabel('Pf concentration (PFU/ml)')
%title('Antibiotic sequestered (%)')



