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

param.MAX_T = 100; % Length of simulation(hour^-1)
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

pMin(16) = 10^(-2); % aMax
pMax(16) = 10^2;

pMin(17) = 10^2; % Starting density for both strains
pMax(17) = 10^8;

n_par = size(pMin,2);

%Log transform min and max
logMin = log(pMin);
logMax = log(pMax);

% Find range spread in log space
diffLog = logMax - logMin;

% Number of data sets generated
reps = 1000;

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


%% 

[sort_pf_rand, sort_no_rand, sort_ratio_rand,p_pf_rand, p_no_rand, p_ratio_rand] = PRCC_function(randomScaled, 1, 0);
[sort_pf_1, sort_no_1, sort_ratio_1, p_pf_1, p_no_1, p_ratio_1] = PRCC_function(randomScaled, 0, 0.1*mic);
[sort_pf_2, sort_no_2, sort_ratio_2, p_pf_2, p_no_2, p_ratio_2] = PRCC_function(randomScaled, 0, 5*mic);
[sort_pf_3, sort_no_3, sort_ratio_3, p_pf_3, p_no_3, p_ratio_3] = PRCC_function(randomScaled, 0, 20*mic);

% Random amax
%%
figure;

subplot(4,3,1)
scatter(1:n_par, table2array(sort_pf_rand))
hold on
xticks([1:n_par])
xticklabels(sort_pf_rand.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel({'Variable AMC','PCC'})
title('Pf+ density')

subplot(4,3,2)
scatter(1:n_par, table2array(sort_no_rand))
hold on
xticks([1:n_par])
xticklabels(sort_no_rand.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel('PCC')
title('Pf- density')

subplot(4,3,3)
scatter(1:n_par, table2array(sort_ratio_rand))
hold on
xticks([1:n_par])
xticklabels(sort_ratio_rand.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel('PCC')
title('Pf+ to Pf- density ratio')

% amax = 0.01

subplot(4,3,4)
scatter(1:n_par, table2array(sort_pf_1))
hold on
xticks([1:n_par])
xticklabels(sort_pf_1.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel({'AMC = 0.1 xMIC','PCC'})
title('Pf+ density')

subplot(4,3,5)
scatter(1:n_par, table2array(sort_no_1))
hold on
xticks([1:n_par])
xticklabels(sort_no_1.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel('PCC')
title('Pf- density')

subplot(4,3,6)
scatter(1:n_par, table2array(sort_ratio_1))
hold on
xticks([1:n_par])
xticklabels(sort_ratio_1.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel('PCC')
title('Pf+ to Pf- density ratio')

% amax = 0.1

subplot(4,3,7)
scatter(1:n_par, table2array(sort_pf_2))
hold on
xticks([1:n_par])
xticklabels(sort_pf_2.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel({'AMC = 5 xMIC','PCC'})
title('Pf+ density')

subplot(4,3,8)
scatter(1:n_par, table2array(sort_no_2))
hold on
xticks([1:n_par])
xticklabels(sort_no_2.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel('PCC')
title('Pf- density')

subplot(4,3,9)
scatter(1:n_par, table2array(sort_ratio_2))
hold on
xticks([1:n_par])
xticklabels(sort_ratio_2.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel('PCC')
title('Pf+ to Pf- density ratio')

% amax = 1

subplot(4,3,10)
scatter(1:n_par, table2array(sort_pf_3))
hold on
xticks([1:n_par])
xticklabels(sort_pf_3.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel({'AMC = 20 xMIC','PCC'})
title('Pf+ density')

subplot(4,3,11)
scatter(1:n_par, table2array(sort_no_3))
hold on
xticks([1:n_par])
xticklabels(sort_no_3.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel('PCC')
title('Pf- density')

subplot(4,3,12)
scatter(1:n_par, table2array(sort_ratio_3))
hold on
xticks([1:n_par])
xticklabels(sort_ratio_3.Properties.RowNames)
xtickangle(0)
ax = gca;
ax.XAxis.FontSize = 8;
ylabel('PCC')
title('Pf+ to Pf- density ratio')



