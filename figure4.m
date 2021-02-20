% Figure 4
% Julie Pourtois
% Last modified: January 8th 2021

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

% 0.1 xMIC
    
param.aMax = 0.1*mic;
n_treatment = total_T/24*nTreatmentPerDay;

param.MAX_T = total_T/n_treatment;

y0_p = [7*10^7 7*10^8 param.aMax];
y0_n = [7*10^7 0 param.aMax];

[t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);

y_p_01_uni = y_p;
t_p_01_uni = t_p;

y_n_01_uni = y_n;
t_n_01_uni = t_n;

% 5 xMIC
    
param.aMax = 5*mic;
n_treatment = total_T/24*nTreatmentPerDay;

param.MAX_T = total_T/n_treatment;

y0_p = [7*10^7 7*10^8 param.aMax];
y0_n = [7*10^7 0 param.aMax];

[t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);

y_p_5_uni = y_p;
t_p_5_uni = t_p;

y_n_5_uni = y_n;
t_n_5_uni = t_n;

% 20 xMIC
    
param.aMax = 20*mic;
n_treatment = total_T/24*nTreatmentPerDay;

param.MAX_T = total_T/n_treatment;

y0_p = [7*10^7 7*10^8 param.aMax];
y0_n = [7*10^7 0 param.aMax];

[t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);

y_p_20_uni = y_p;
t_p_20_uni = t_p;

y_n_20_uni = y_n;
t_n_20_uni = t_n;

% 30 xMIC
    
param.aMax = 30*mic;
n_treatment = total_T/24*nTreatmentPerDay;

param.MAX_T = total_T/n_treatment;

y0_p = [7*10^7 7*10^8 param.aMax];
y0_n = [7*10^7 0 param.aMax];

[t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);

y_p_30_uni = y_p;
t_p_30_uni = t_p;

y_n_30_uni = y_n;
t_n_30_uni = t_n;

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

threshold = 1;

equiDens_005 = zeros(reps,2);
equiDens_01 = zeros(reps,2);
equiDens_1 = zeros(reps,2);
equiDens_2 = zeros(reps,2);
equiDens_5 = zeros(reps,2);
equiDens_8 = zeros(reps,2);
equiDens_20 = zeros(reps,2);
equiDens_30 = zeros(reps,2);

over_count_005 = 0;
over_count_01 = 0;
over_count_1 = 0;
over_count_2 = 0;
over_count_5 = 0;
over_count_8 = 0;
over_count_20 = 0;
over_count_30 = 0;

extinct_005 = 0;
extinct_01 = 0;
extinct_1 = 0;
extinct_2 = 0;
extinct_5 = 0;
extinct_8 = 0;
extinct_20 = 0;
extinct_30 = 0;

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
    
    % 0.05 xMIC
   
    param.aMax = 0.05*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 | y_n_last_day > 1;
    
    equiDens_005(i,1) = y_p_last_day;
    equiDens_005(i,2) = y_n_last_day;
    
    if (y_n_last_day/y_p_last_day > threshold & no_extinct)
        
        over_count_005 = over_count_005 + 1;
        over_005(over_count_005,1) = y_p_last_day;
        over_005(over_count_005,2) = y_n_last_day;
        
        over_param_005(over_count_005,:) = randomScaled(i,:);
        
        
    end
    
    if (no_extinct == 0)        
        extinct_005 = extinct_005 + 1;
    end
    
    % 0.1 xMIC
   
    param.aMax = 0.1*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_01{i} = y_p;
    t_p_01{i} = t_p;
    
    y_n_01{i} = y_n;
    t_n_01{i} = t_n;
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
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
    
    % 1 xMIC
    
    param.aMax = mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    equiDens_1(i,1) = y_p_last_day;
    equiDens_1(i,2) = y_n_last_day;
    
    no_extinct = (y_p_last_day > 1 | y_n_last_day > 1);
    
    if (y_n_last_day/y_p_last_day > threshold & no_extinct)
        
        over_count_1 = over_count_1 + 1;
        over_1(over_count_1,1) = y_p_last_day;
        over_1(over_count_1,2) = y_n_last_day;
        
        over_param_1(over_count_1,:) = randomScaled(i,:);
        
        
    end
    
    if (no_extinct == 0)        
        extinct_1 = extinct_1 + 1;
    end
    
    % 2 xMIC
  
    param.aMax = 2*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    equiDens_2(i,1) = y_p_last_day;
    equiDens_2(i,2) = y_n_last_day;
    
    no_extinct = y_p_last_day > 1 | y_n_last_day > 1;
    
    if (y_n_last_day/y_p_last_day > threshold & no_extinct)
        
        over_count_2 = over_count_2 + 1;
        over_2(over_count_2,1) = y_p_last_day;
        over_2(over_count_2,2) = y_n_last_day;
        
        over_param_2(over_count_2,:) = randomScaled(i,:);
        
        
    end
    
    if (no_extinct == 0)        
        extinct_2 = extinct_2 + 1;
    end
    
    % 5 xMIC
    
    param.aMax = 5*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_5{i} = y_p;
    t_p_5{i} = t_p;
    
    y_n_5{i} = y_n;
    t_n_5{i} = t_n;
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
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
    
    % 8 xMIC
    
    param.aMax = 8*mic;
    
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    equiDens_8(i,1) = y_p_last_day;
    equiDens_8(i,2) = y_n_last_day;
    
    no_extinct = y_p_last_day > 1 | y_n_last_day > 1;
    
    if (y_n_last_day/y_p_last_day > threshold & no_extinct)
        
        over_count_8 = over_count_8 + 1;
        over_8(over_count_8,1) = y_p_last_day;
        over_8(over_count_8,2) = y_n_last_day;
        
        over_param_8(over_count_8,:) = randomScaled(i,:);
        
        
    end
    
    if (no_extinct == 0)        
        extinct_8 = extinct_8 + 1;
    end
    
    % 20 xMIC
    
    param.aMax = 20*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_20{i} = y_p;
    t_p_20{i} = t_p;
    
    y_n_20{i} = y_n;
    t_n_20{i} = t_n;
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
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
    
    
    % 30 xMIC
    
    param.aMax = 30*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_30{i} = y_p;
    t_p_30{i} = t_p;
    
    y_n_30{i} = y_n;
    t_n_30{i} = t_n;
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    equiDens_30(i,1) = y_p_last_day;
    equiDens_30(i,2) = y_n_last_day;
    
    no_extinct = y_p_last_day > 1 | y_n_last_day > 1;
    
    if (y_n_last_day/y_p_last_day > threshold & no_extinct)
        
        
        over_count_30 = over_count_30 + 1;
        
        over_30(over_count_30,1) = y_p_last_day;
        over_30(over_count_30,2) = y_n_last_day;
        
        over_param_30(over_count_30,:) = randomScaled(i,:);
        
        
    end
    
    if (no_extinct == 0)        
        extinct_30 = extinct_30 + 1;
    end
    
    
end

%% Convert from reps to percent

to_percent = reps/100;

under_count_005 = (reps - over_count_005 - extinct_005)/to_percent;
under_count_01 = (reps - over_count_01 - extinct_01)/to_percent;
under_count_1 = (reps - over_count_1 - extinct_1)/to_percent;
under_count_2 = (reps - over_count_2 - extinct_2)/to_percent;
under_count_5 = (reps - over_count_5 - extinct_5)/to_percent;
under_count_8 = (reps - over_count_8 - extinct_8)/to_percent;
under_count_20 = (reps - over_count_20 - extinct_20)/to_percent;
under_count_30 = (reps - over_count_30 - extinct_30)/to_percent;

over_count_01 = over_count_01/to_percent;
over_count_5 = over_count_5/to_percent;
over_count_20 = over_count_20/to_percent;
over_count_30 = over_count_30/to_percent;

extinct_01 = extinct_01/to_percent;
extinct_5 = extinct_5/to_percent;
extinct_20 = extinct_20/to_percent;
extinct_30 = extinct_30/to_percent;

%%
figure()
scatter(equiDens_30(:,2),equiDens_30(:,1), 60, 'filled', 'MarkerFaceAlpha',.3)
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [1, 10^10];
y = [1, 10^10];
%xticks([10^7 10^8])
daspect([1 1 1])
line(x,y,'Color','red','LineWidth',1.5)
text(2,10^9,sprintf('%.0f / %.0f / %.0f', under_count_30,over_count_30,extinct_30))%[under_count_005],[over_count_005],[extinct_005])
xlabel('Pf- density (CFU/ml)')
ylabel('Pf+ density (CFU/ml)')
title('30 xMIC')

figure()
plot(t_p_30{1},y_p_30{1}(:,1),'Color',[0.8, 0.8, 0.8])
hold on
ylim([1 5*10^9])
xlim([1 120])
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Pf+ density (CFU/ml)')
title('30 xMIC')

for j = 1:size(t_p_30,2)
      
    plot(t_p_30{j},y_p_30{j}(:,1),'Color',[0.8, 0.8, 0.8,0.4])  
    
end

plot(t_p_30_uni,y_p_30_uni(:,1), 'LineWidth', 2, 'Color', [0 0.4470 0.7410])


%% 
figure()

subplot(3,3,1)
scatter(equiDens_01(:,2),equiDens_01(:,1), 60, 'filled', 'MarkerFaceAlpha',.3)
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [5*10^6, 2*10^8];
y = [5*10^6, 2*10^8];
xticks([10^7 10^8])
daspect([1 1 1])
line(x,y,'Color','red','LineWidth',1.5)
text(7*10^6,1.5*10^8,sprintf('%.0f / %.0f / %.0f', under_count_01,over_count_01,extinct_01))%[under_count_005],[over_count_005],[extinct_005])
xlabel('Pf- density (CFU/ml)')
ylabel('Pf+ density (CFU/ml)')
title('0.1 xMIC')

subplot(3,3,2)
plot(t_p_01{1},y_p_01{1}(:,1),'Color',[0.8, 0.8, 0.8])
hold on
ylim([1 5*10^8])
xlim([1 120])
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Pf+ density (CFU/ml)')
title('0.1 xMIC')

for j = 1:size(t_p_01,2)
      
    plot(t_p_01{j},y_p_01{j}(:,1),'Color',[0.8, 0.8, 0.8,0.4])  
    
end

plot(t_p_01_uni,y_p_01_uni(:,1), 'LineWidth', 2, 'Color', [0 0.4470 0.7410])


subplot(3,3,3)
plot(t_n_01{1},y_n_01{1}(:,1),'Color',[0.8, 0.8, 0.8])
hold on
ylim([1 5*10^8])
xlim([1 120])
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Pf- density (CFU/ml)')
title('0.1 xMIC')

for j = 1:size(t_n_01,2)
      
    plot(t_n_01{j},y_n_01{j}(:,1),'Color',[0.8, 0.8, 0.8])  
    
end

plot(t_n_01_uni,y_n_01_uni(:,1), 'LineWidth', 2, 'Color', [0 0.4470 0.7410])


subplot(3,3,4)
scatter(equiDens_20(:,2),equiDens_20(:,1), 60, 'filled','MarkerFaceAlpha',.3)
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [1, 10^10];
y = [1, 10^10];
daspect([1 1 1])
line(x,y,'Color','red','LineWidth',1.5)
text(2,10^9,sprintf('%.0f / %.0f / %.0f', under_count_20,over_count_20,extinct_20))%[under_count_005],[over_count_005],[extinct_005])
xlabel('Pf- density (CFU/ml)')
ylabel('Pf+ density (CFU/ml)')
title('20 xMIC')

subplot(3,3,5)
plot(t_p_20{1},y_p_20{1}(:,1),'Color',[0.8, 0.8, 0.8])
hold on
ylim([1 5*10^8])
xlim([1 120])
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Pf+ density (CFU/ml)')
title('20 xMIC')

for j = 1:size(t_p_20,2)
      
    plot(t_p_20{j},y_p_20{j}(:,1),'Color',[0.8, 0.8, 0.8,0.4])  
    
end

plot(t_p_20_uni,y_p_20_uni(:,1), 'LineWidth', 2, 'Color', [0 0.4470 0.7410])

subplot(3,3,6)
plot(t_n_20{1},y_n_20{1}(:,1),'Color',[0.8, 0.8, 0.8])
hold on
ylim([1 5*10^8])
xlim([1 120])
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Pf- density (CFU/ml)')
title('20 xMIC')

for j = 1:size(t_n_20,2)
      
    plot(t_n_20{j},y_n_20{j}(:,1),'Color',[0.8, 0.8, 0.8,0.4])  
    
end

plot(t_n_20_uni,y_n_20_uni(:,1), 'LineWidth', 2, 'Color', [0 0.4470 0.7410])



subplot(3,3,7)
s = scatter(equiDens_30(:,2),equiDens_30(:,1), 60, 'filled','MarkerFaceAlpha',.3)
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [1, 10^10];
y = [1, 10^10];
daspect([1 1 1])
line(x,y,'Color','red','LineWidth',1.5)
text(2,10^9,sprintf('%.0f / %.0f / %.0f', under_count_30,over_count_30,extinct_30))%[under_count_005],[over_count_005],[extinct_005])
xlabel('Pf- density (CFU/ml)')
ylabel('Pf+ density (CFU/ml)')
title('30 xMIC')


%t_val = (1:240)./2;


subplot(3,3,8)
plot(t_p_30{1},y_p_30{1}(:,1),'Color',[0.8, 0.8, 0.8])
hold on
ylim([1 5*10^9])
xlim([1 120])
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Pf+ density (CFU/ml)')
title('30 xMIC')

for j = 1:size(t_p_30,2)
      
    plot(t_p_30{j},y_p_30{j}(:,1),'Color',[0.8, 0.8, 0.8,0.4])  
    
end

plot(t_p_30_uni,y_p_30_uni(:,1), 'LineWidth', 2, 'Color', [0 0.4470 0.7410])



%t_val = (1:240)./2;


subplot(3,3,9)
plot(t_n_30{1},y_n_30{1}(:,1),'Color',[0.8, 0.8, 0.8])
hold on
ylim([1 5*10^9])
xlim([1 120])
set(gca,'yscale','log')
xlabel('Time (hours)')
ylabel('Pf- density (CFU/ml)')
title('30 xMIC')

for j = 1:size(t_n_30,2)
      
    plot(t_n_30{j},y_n_30{j}(:,1),'Color',[0.8, 0.8, 0.8,0.4])  
    
end

plot(t_n_30_uni,y_n_30_uni(:,1), 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
