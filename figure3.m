% Figure 3
% General sensitivity analysis, direct competition
% Last modified: August 19th 2020, Julie Pourtois
%%

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
param.phi = 10^6; % Antibiotic sequestration factor

mic = param.kc/(-param.rMin/(param.rMax - param.deltaB) - 1)^(1/param.hill);

total_T = 120; % Total duration of treatment (hour)
nTreatmentPerDay = 2; % Number treatments per day

%% General sensitivity

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

pMin(16) = 10^2; % Starting density for both strains
pMax(16) = 10^8;

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

threshold = 1;

equiDens_01 = zeros(reps,2);
equiDens_5 = zeros(reps,2);
equiDens_20 = zeros(reps,2);

over_count_01 = 0;
over_count_5 = 0;
over_count_20 = 0;

extinct_01 = 0;
extinct_5 = 0;
extinct_20 = 0;

for i = 1:reps
    
    param.hill = randomScaled(i,1); % Hill parameter
    param.rMax = randomScaled(i,2); % Growth rate (hour^-1)
    param.K = 10^8; % Carrying capacity
    param.rMin = -randomScaled(i,3); % Maximum killing rate
    param.epsiR = randomScaled(i,4); % Growth rate depedence of kill rate
    param.gammaR = randomScaled(i,5); % Growth rate for half kill rate % 0.1
    param.w = 468*10^6; % Antibiotic molecular weight
    param.epsiK = randomScaled(i,6); % Density-dependence of kill rate
    param.kc = 3.68; % k0 in manuscript. Antibiotic concentration at which killing is half of max
    param.xi = randomScaled(i,7); % Max increase density-dependent kill
    param.gammaK = randomScaled(i,8); % Density for k at half of its maximum
    param.deltaB = randomScaled(i,9); % Basal death rate (hour^-1)
    param.deltaA = randomScaled(i,10); % Antibiotic decay rate (hour^-1) % 0.05
    param.lambda = randomScaled(i,11); %Phage production rate (hour^-1)
    param.deltaV = randomScaled(i,12); %Phage decay rate (hour^-1)
    param.kD = randomScaled(i,13);% Binding dissociation factor 10^13-10^15
    param.exponent = 1;
    param.theta = randomScaled(i,14); % Metabolic cost
    param.pDeath = 1; % Additional mortality from phage production
    param.phi = randomScaled(i,15); % Antibiotic sequestration factor
    param.init = randomScaled(i,16); % Starting concentration of both strains
    
    total_T = 120; % Total duration of treatment (hour)
    nTreatmentPerDay = 2; % Number treatments per day
    
    
    % 0.1 xMIC
   
    param.aMax = 0.1*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0 = [7*10^7 7*10^7 7*10^8 param.aMax];
    
    [t,y] = competeRegimens(n_treatment, y0, param);
    
    y_01{i} = y;
    t_01{i} = t;
    
    y_p_last_day = median(y(t > 96,1));
    y_n_last_day = median(y(t > 96,2));
    
    equiDens_01(i,1) = y_p_last_day;
    equiDens_01(i,2) = y_n_last_day;
    
    no_extinct = y_p_last_day > 1 | y_n_last_day > 1;
    
    if (y_n_last_day/y_p_last_day > threshold & no_extinct)
        
        over_count_01 = over_count_01 + 1;
        over_01(over_count_01,1) = y_p_last_day;
        over_01(over_count_01,2) = y_n_last_day;
        
        over_param_01(over_count_01,:) = randomScaled(i,:);
        
        
    end
    
    if (no_extinct == 0)        
        extinct_01 = extinct_01 + 1;
    end
    
    % 5 xMIC
    
    param.aMax = 5*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0 = [7*10^7 7*10^7 7*10^8 param.aMax];
    
    [t,y] = competeRegimens(n_treatment, y0, param);
    
    y_5{i} = y;
    t_5{i} = t;
    
    y_p_last_day = median(y(t > 96,1));
    y_n_last_day = median(y(t > 96,2));
    
    equiDens_5(i,1) = y_p_last_day;
    equiDens_5(i,2) = y_n_last_day;
    
    no_extinct = y_p_last_day > 1 | y_n_last_day > 1;
    
    if (y_n_last_day/y_p_last_day > threshold & no_extinct)
        
        over_count_5 = over_count_5 + 1;
        over_5(over_count_5,1) = y_p_last_day;
        over_5(over_count_5,2) = y_n_last_day;
        
        over_param_5(over_count_5,:) = randomScaled(i,:);
        
        
    end
    
    if (no_extinct == 0)        
        extinct_5 = extinct_5 + 1;
    end
    
    % 20 xMIC
    
    param.aMax = 20*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0 = [7*10^7 7*10^7 7*10^8 param.aMax];
    
    [t,y] = competeRegimens(n_treatment, y0, param);
    
    y_20{i} = y;
    t_20{i} = t;
    
    y_p_last_day = median(y(t > 96,1));
    y_n_last_day = median(y(t > 96,2));
    
    equiDens_20(i,1) = y_p_last_day;
    equiDens_20(i,2) = y_n_last_day;
    
    no_extinct = y_p_last_day > 1 | y_n_last_day > 1;
    
    if (y_n_last_day/y_p_last_day > threshold & no_extinct)
        
        
        over_count_20 = over_count_20 + 1;
        
        over_20(over_count_20,1) = y_p_last_day;
        over_20(over_count_20,2) = y_n_last_day;
        
        over_param_20(over_count_20,:) = randomScaled(i,:);
        
        
    end
    
    if (no_extinct == 0)        
        extinct_20 = extinct_20 + 1;
    end
end


%%

under_count_01 = reps - over_count_01 - extinct_01;
under_count_5 = reps - over_count_5 - extinct_5;
under_count_20 = reps - over_count_20 - extinct_20;

%%

figure()

subplot(2,3,4)
scatter(equiDens_01(:,2),equiDens_01(:,1), 60, 'filled', 'MarkerFaceAlpha',.3)
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [10^4 10^9];
y = [10^4 10^9];
daspect([1 1 1])
xlim([10^4 10^9])
ylim([10^4 10^9])
line(x,y,'Color','red','LineWidth',1.5)
text(2*10^4,3*10^8,sprintf('%.0f / %.0f / %.0f', under_count_01,over_count_01,extinct_01))%[under_count_005],[over_count_005],[extinct_005])
xlabel('Pf- density (CFU/ml)')
ylabel('Pf+ density (CFU/ml)')
title('0.1 xMIC')

subplot(2,3,5)
scatter(equiDens_5(:,2),equiDens_5(:,1), 60, 'filled','MarkerFaceAlpha',.3)
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [1, 10^10];
y = [1, 10^10];
daspect([1 1 1])
line(x,y,'Color','red','LineWidth',1.5)
text(2,10^9,sprintf('%.0f / %.0f / %.0f', under_count_5,over_count_5,extinct_5))%[under_count_005],[over_count_005],[extinct_005])
xlabel('Pf- density (CFU/ml)')
title('5 xMIC')

subplot(2,3,6)
s = scatter(equiDens_20(:,2),equiDens_20(:,1), 60, 'filled','MarkerFaceAlpha',.3)
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [1, 10^10];
y = [1, 10^10];
daspect([1 1 1])
line(x,y,'Color','red','LineWidth',1.5)
text(2,10^9,sprintf('%.0f / %.0f / %.0f', under_count_20,over_count_20,extinct_20))%[under_count_005],[over_count_005],[extinct_005])
xlabel('Pf- density (CFU/ml)')
title('20 xMIC')


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


%%% Figure 3A: Death rate versus AMC - No dependence on growth rate 

AC = round(10.^(-2:0.25:2.5),2);
B= 5*10^6;
V = 5*10^7;

r_p = zeros(length(AC),1); % Pf+
r_n = zeros(length(AC),1); % Pf-

param3A = param;
param3A.epsiR = 0;
param3A.theta = 0.5;

for i = 1:length(AC)
    
    r_p(i) = growthRate(B,V,AC(i),param3A,1);
    r_n(i) = growthRate(B,V,AC(i),param3A,0);
end


subplot(2,3,1)
plot(AC/mic, -r_p)
hold on
plot(AC/mic, -r_n)
legend('Pf (+)', 'Pf (-)','Location','northwest')
set(gca,'xscale','log')
xlabel('Antibiotic concentration (xMIC)')
ylabel('Death rate (h^{-1})')
xlim([10^(-1) 10^(3)])

%%% Figure 3B: Death rate versus AMC - Dependence on growth rate 

param3B = param;
param3B.epsiR = 1;
param3B.theta = 0.5;

B= 5*10^7;
V = 5*10^8;

% Phage-positive
r_p = zeros(length(AC),1);
r_n = zeros(length(AC),1);

for i = 1:length(AC)
    
    r_p(i) = growthRate(B,V,AC(i),param3B,1);
    r_n(i) = growthRate(B,V,AC(i),param3B,0);
end

subplot(2,3,2)
plot(AC/mic, -r_p)
hold on
plot(AC/mic, -r_n)
legend('Pf (+)', 'Pf (-)','Location','northwest')
set(gca,'xscale','log')
xlabel('Antibiotic concentration (xMIC)')
ylabel('Death rate (h^{-1})')
xlim([10^(-1) 10^(3)])

%%% Figure 3C - Dependence of growth rate on metabolic rate 

r_diff = zeros(length(AC),3);
thetas = [0.2,0.5,0.8];

B= 5*10^7;
V = 5*10^8;

param3C = param;

for l = 1:3
    
    param3C.theta = thetas(l);
    
    r_p = zeros(length(AC),1); % Pf+
    r_n = zeros(length(AC),1); % Pf-
    
    for i = 1:length(AC)
        
        r_p(i) = growthRate(B,V,AC(i),param3C,1);
        r_n(i) = growthRate(B,V,AC(i),param3C,0);
    end
    
    r_diff(:,l) = r_p - r_n;
    
end

subplot(2,3,3)
plot(AC/mic, -r_diff(:,1))
hold on
plot(AC/mic, -r_diff(:,2))
plot(AC/mic, -r_diff(:,3))
set(gca,'xscale','log')
xlabel('Antibiotic concentration (xMIC)')
ylabel('Effect of Pf on death rate (h^{-1})')
legend('\theta = 0.2', '\theta = 0.5', '\theta = 0.8', 'Location','southwest')
xlim([10^(-1) 10^(3)])
