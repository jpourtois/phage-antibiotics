% Figure 1B,C,D
% Julie Pourtois
% Last modified: Jan 7th 2021

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

% Calculate MIC

mic = param.kc/(-param.rMin/(param.rMax - param.deltaB) - 1)^(1/param.hill);

% Default value for varying parameters 

param.theta = 0.2; % Metabolic cost
param.phi = 10^6; % Antibiotic sequestration factor
param.aMax = 3*mic; % Peak antibiotic concentration



%%

total_T = 120;
nTreatmentPerDay = 2; % Number treatments per day

n_treatment = total_T/24*nTreatmentPerDay;

param.MAX_T = total_T/n_treatment;

y0_p = [10^7 10^8 param.aMax];
y0_n = [10^7 0 param.aMax];

[t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);

%%

figure;
subplot(1,3,1)
plot(t_p, y_p(:,1),'Color', [0.8, 0.4, 0.4], 'LineWidth', 2)
set(gca, 'YScale', 'log')
ylim([5*10^6 10^8])
xlim([0 120])
xlabel('Time (hours)')
ylabel('Pf+ bacterial density (CFU/ml)')


subplot(1,3,2)
plot(t_p, y_p(:,2),'Color', [0.8, 0.4, 0.4], 'LineWidth', 2)
set(gca, 'YScale', 'log')
ylim([5*10^7 10^9])
xlim([0 120])
xlabel('Time (hours)')
ylabel('Pf density (PFU/ml)')

subplot(1,3,3)
plot(t_p, y_p(:,3)/mic,'Color', [0.8, 0.4, 0.4], 'LineWidth', 2)
ylim([0 4])
xlim([0 120])
xlabel('Time (hours)')
ylabel('Antibiotic concentration (xMIC)')



