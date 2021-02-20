%%% 
% Figure 5
% Last modified: January 11th 2021, Julie Pourtois
%

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


%%

[percent_win_theta,percent_win_n_theta, range2_theta] = fc_var_par_median(randomScaled, 'theta',1);
[percent_win_phi, percent_win_n_phi, range2_phi] = fc_var_par_median(randomScaled, 'phi',1);
[percent_win_deltaV, percent_win_n_deltaV, range2_deltaV] = fc_var_par_median(randomScaled, 'deltaV',1);
[percent_win_deltaA, percent_win_n_deltaA, range2_deltaA] = fc_var_par_median(randomScaled, 'deltaA',1);
[percent_win_rMin,percent_win_n_rMin, range2_rMin] = fc_var_par_median(randomScaled, 'rMin',1);
[percent_win_epsiR,percent_win_n_epsiR, range2_epsiR] = fc_var_par_median(randomScaled, 'epsiR',1);

%%
range_dose = 10.^(-2:0.05:2);

figure;
subplot(2,3,1)
s = surface(range2_theta, range_dose./mic, percent_win_theta, 'FaceAlpha',1)
colorbar
hold on
%shading interp
s.EdgeColor = 'none';
set(gca,'yscale','log')
ylim([10^(-1) 10^(2.8)])
xlim([0.1 0.9])
xlabel('Metabolic cost \theta')
ylabel('AMC (xMIC)')

subplot(2,3,2)
s = surface(range2_phi, range_dose./mic, percent_win_phi, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
set(gca,'xscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlabel('Antibiotic sequestration \phi')
ylabel('AMC (xMIC)')
%title('Proportion of runs won by Pf+')


subplot(2,3,3)
s = surface(range2_deltaV, range_dose./mic, percent_win_deltaV, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
set(gca,'xscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlabel('Phade decay \delta_V')
ylabel('AMC (xMIC)')

subplot(2,3,4)
s = surface(range2_epsiR, range_dose./mic, percent_win_epsiR, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlim([0.1 0.9])
xlabel("Replication-dependence \epsilon_R")
ylabel('AMC (xMIC)')

subplot(2,3,5)
s = surface(range2_deltaA, range_dose./mic, percent_win_deltaA, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlim([0.05 0.5])
xlabel('Antibiotic decay rate \delta_A')
ylabel('AMC (xMIC)')

subplot(2,3,6)
s = surface(range2_rMin, range_dose./mic, percent_win_rMin, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlim([2 15])
xlabel('Antibiotic max killing rate \Gamma')
ylabel('AMC (xMIC)')


%%

range_dose = 10.^(-2:0.05:2);

figure;
subplot(2,3,1)
s = surface(range2_theta, range_dose./mic, percent_win_n_theta, 'FaceAlpha',1)
colorbar
hold on
%shading interp
s.EdgeColor = 'none';
set(gca,'yscale','log')
ylim([10^(-1) 10^(2.8)])
xlim([0.1 0.9])
xlabel('Metabolic cost \theta')
ylabel('AMC (xMIC)')

subplot(2,3,2)
s = surface(range2_phi, range_dose./mic, percent_win_n_phi, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
set(gca,'xscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlabel('Antibiotic sequestration \phi')
ylabel('AMC (xMIC)')
%title('Proportion of runs won by Pf+')


subplot(2,3,3)
s = surface(range2_deltaV, range_dose./mic, percent_win_n_deltaV, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
set(gca,'xscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlabel('Phade decay \delta_V')
ylabel('AMC (xMIC)')

subplot(2,3,4)
s = surface(range2_epsiR, range_dose./mic, percent_win_n_epsiR, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlim([0.1 0.9])
xlabel("Replication-dependence \epsilon_R")
ylabel('AMC (xMIC)')

subplot(2,3,5)
s = surface(range2_deltaA, range_dose./mic, percent_win_n_deltaA, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlim([0.05 0.5])
xlabel('Antibiotic decay rate \delta_A')
ylabel('AMC (xMIC)')

subplot(2,3,6)
s = surface(range2_rMin, range_dose./mic, percent_win_n_rMin, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlim([2 15])
xlabel('Antibiotic max killing rate \Gamma')
ylabel('AMC (xMIC)')



%% 
range_dose = 10.^(-2:0.05:2);
figure;
subplot(2,3,1)
s = surface(range2_theta, range_dose./mic, percent_win_theta - percent_win_n_theta, 'FaceAlpha',1)
colorbar
hold on
%shading interp
s.EdgeColor = 'none';
set(gca,'yscale','log')
ylim([10^(-1) 10^(2.8)])
xlim([0.1 0.9])
xlabel('Metabolic cost \theta')
ylabel('AMC (xMIC)')

subplot(2,3,2)
s = surface(range2_phi, range_dose./mic, percent_win_phi - percent_win_n_phi, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
set(gca,'xscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlabel('Antibiotic sequestration \phi')
ylabel('AMC (xMIC)')
%title('Proportion of runs won by Pf+')


subplot(2,3,3)
s = surface(range2_deltaV, range_dose./mic, percent_win_deltaV - percent_win_n_deltaV, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
set(gca,'xscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlabel('Phade decay \delta_V')
ylabel('AMC (xMIC)')

subplot(2,3,4)
s = surface(range2_deltaA, range_dose./mic, percent_win_deltaA - percent_win_n_deltaA, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlim([0.05 0.5])
xlabel('Antibiotic decay rate \delta_A')
ylabel('AMC (xMIC)')

subplot(2,3,5)
s = surface(range2_rMin, range_dose./mic, percent_win_rMin - percent_win_n_rMin, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlim([2 15])
xlabel('Antibiotic max killing rate \Gamma')
ylabel('AMC (xMIC)')

subplot(2,3,6)
s = surface(range2_epsiR, range_dose./mic, percent_win_epsiR - percent_win_n_epsiR, 'FaceAlpha',1)
colorbar
hold on
%shading interp
set(gca,'yscale','log')
s.EdgeColor = 'none';
ylim([10^(-1) 10^(2.8)])
xlim([0.1 0.9])
xlabel("Replication-dependence \epsilon_R")
ylabel('AMC (xMIC)')



%% Export for analysis in R

diff_theta = percent_win_theta - percent_win_n_theta;
diff_phi = percent_win_phi - percent_win_n_phi;
diff_deltaV = percent_win_deltaV - percent_win_n_deltaV;
diff_epsiR = percent_win_epsiR - percent_win_n_epsiR;
diff_deltaA = percent_win_deltaA - percent_win_n_deltaA;
diff_rMin = percent_win_rMin - percent_win_n_rMin;


%%

writematrix(diff_theta,'diff_theta.csv')
writematrix(diff_phi,'diff_phi.csv')
writematrix(diff_deltaV, 'diff_deltaV.csv')
writematrix(diff_epsiR,'diff_epsiR.csv')
writematrix(diff_deltaA, 'diff_deltaA.csv')
writematrix(diff_rMin, 'diff_rMin.csv')






