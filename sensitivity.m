%%% General sensitivity analysis
% Generate supplemental figures S1 and S2
% Last modified: August 19th 2020, Julie Pourtois
%
%% General sensitivity

pMin(1) = 0.6; % Hill 
pMax(1) = 2.5;

pMin(2) = 0.5; % rMax
pMax(2) = 1.5;

pMin(3) = 2; % rMin
pMax(3) = 15;

pMin(4)= 0.05; % gammaR
pMax(4) = 0.5;

pMin(5) = 10^(-1); % epsiK
pMax(5) = 1;

pMin(6) = 2; % xi
pMax(6) = 20;

pMin(7) = 10^6; % gammaK
pMax(7) = 10^7;

pMin(8) = 0.025; % deltaB
pMax(8) = 0.25;


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
randomUnscaled = lhsdesign(reps,8);

% Scale the random number with the range spread in log space
randomdiffLog = zeros(reps,8);
 
for n = 1:8
    randomdiffLog(:,n) = randomUnscaled(:,n)*diffLog(n);    
end

randomLogScaled = logMin + randomdiffLog;

%Get out of log space
randomScaled = exp(randomLogScaled);

threshold = 1;

equiDens_005 = zeros(reps,2);
equiDens_05 = zeros(reps,2);
equiDens_1 = zeros(reps,2);
equiDens_2 = zeros(reps,2);
equiDens_5 = zeros(reps,2);
equiDens_8 = zeros(reps,2);
equiDens_15 = zeros(reps,2);

over_count_005 = 0;
over_count_05 = 0;
over_count_1 = 0;
over_count_2 = 0;
over_count_5 = 0;
over_count_8 = 0;
over_count_15 = 0;

extinct_005 = 0;
extinct_05 = 0;
extinct_1 = 0;
extinct_2 = 0;
extinct_5 = 0;
extinct_8 = 0;
extinct_15 = 0;

for i = 1:reps
    
    param.hill = randomScaled(i,1); % Hill parameter
    param.rMax = randomScaled(i,2); % Growth rate (hour^-1)
    param.K = 10^8; % Carrying capacity
    param.rMin = -randomScaled(i,3); % Maximum killing rate
    param.epsiR = 0.5; % Growth rate depedence of kill rate
    param.gammaR = randomScaled(i,4); % Growth rate for half kill rate % 0.1
    param.w = 468*10^6; % Antibiotic molecular weight
    param.epsiK = randomScaled(i,5); % Density-dependence of kill rate
    param.kc = 3.68; % k0 in manuscript. Antibiotic concentration at which killing is half of max
    param.xi = randomScaled(i,6); % Max increase density-dependent kill
    param.gammaK = randomScaled(i,7); % Density for k at half of its maximum
    param.deltaB = randomScaled(i,8); % Basal death rate (hour^-1)
    param.deltaA = 0.1; % Antibiotic decay rate (hour^-1) % 0.05
    param.lambda = 1; %Phage production rate (hour^-1)
    param.deltaV = 0.1; %Phage decay rate (hour^-1)
    param.kD = 10^(14);% Binding dissociation factor 10^13-10^15
    param.exponent = 1;
    param.theta = 0.2; % Metabolic cost
    param.pDeath = 1; % Additional mortality from phage production
    param.phi = 10^6; % Antibiotic sequestration factor
    param.aMax = 1; % Peak antibiotic concentration
      
    total_T = 120; % Total duration of treatment (hour)
    nTreatmentPerDay = 1; % Number treatments per day
    
    % Calculate MIC

    mic = param.kc/(-param.rMin/(param.rMax - param.deltaB) - 1)^(1/param.hill);
    
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
    
    % 0.5 xMIC
   
    param.aMax = 0.5*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    equiDens_05(i,1) = y_p_last_day;
    equiDens_05(i,2) = y_n_last_day;
    
    no_extinct = y_p_last_day > 1 | y_n_last_day > 1;
    
    if (y_n_last_day/y_p_last_day > threshold & no_extinct)
        
        over_count_05 = over_count_05 + 1;
        over_05(over_count_05,1) = y_p_last_day;
        over_05(over_count_05,2) = y_n_last_day;
        
        over_param_05(over_count_05,:) = randomScaled(i,:);
        
        
    end
    
    if (no_extinct == 0)        
        extinct_05 = extinct_05 + 1;
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
    aPerDay = 8*mic ; % Total antibiotics dose per day % 1
    
    param.aMax = aPerDay/nTreatmentPerDay;
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
    
    % 15 xMIC
    
    param.aMax = 15*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    equiDens_15(i,1) = y_p_last_day;
    equiDens_15(i,2) = y_n_last_day;
    
    no_extinct = y_p_last_day > 1 | y_n_last_day > 1;
    
    if (y_n_last_day/y_p_last_day > threshold & no_extinct)
        
        
        over_count_15 = over_count_15 + 1;
        
        over_15(over_count_15,1) = y_p_last_day;
        over_15(over_count_15,2) = y_n_last_day;
        
        over_param_15(over_count_15,:) = randomScaled(i,:);
        
        
    end
    
    if (no_extinct == 0)        
        extinct_15 = extinct_15 + 1;
    end
end

%%

range = log10(pMax) - log10(pMin);

for i = 1:length(over_param_2(:,1))
    
    min_param = log10(over_param_2(i,:)) - log10(pMin);
    
    dist(i,:) = min_param./range;
    
end

figure();
boxplot(dist)
yticks([0 1])
yticklabels({'Min','Max'})
format_ticks(gca,{'H','r_{max}','\Gamma','\gamma_{r}','\epsilon_{k}','\xi','\gamma _{k}','\delta _b'},[],[],[],[],[],[],[],'FontSize',12)

%%

under_count_005 = reps - over_count_005 - extinct_005;
under_count_05 = reps - over_count_05 - extinct_05;
under_count_1 = reps - over_count_1 - extinct_1;
under_count_2 = reps - over_count_2 - extinct_2;
under_count_5 = reps - over_count_5 - extinct_5;
under_count_8 = reps - over_count_8 - extinct_8;
under_count_15 = reps - over_count_15 - extinct_15;

%% 
figure()
subplot(2,3,1)
scatter(equiDens_005(:,2),equiDens_005(:,1))
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [10^7, 2*10^8];
y = [10^7, 2*10^8];
line(x,y,'Color','red','LineWidth',1.5)
text(1.2*10^7,1.5*10^8,sprintf('%.0f / %.0f / %.0f', under_count_005,over_count_005,extinct_005))%[under_count_005],[over_count_005],[extinct_005])
ylabel('Pf+ density (CFU/ml)')
title('0.05 xMIC')

subplot(2,3,2)
scatter(equiDens_05(:,2),equiDens_05(:,1))
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [10^7, 2*10^8];
y = [10^7, 2*10^8];
line(x,y,'Color','red','LineWidth',1.5)
text(1.2*10^7,1.5*10^8,sprintf('%.0f / %.0f / %.0f', under_count_05,over_count_05,extinct_05))%[under_count_005],[over_count_005],[extinct_005])
title('0.5 xMIC')

subplot(2,3,3)
scatter(equiDens_1(:,2),equiDens_1(:,1))
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [10^7, 2*10^8];
y = [10^7, 2*10^8];
line(x,y,'Color','red','LineWidth',1.5)
text(1.2*10^7,1.5*10^8,sprintf('%.0f / %.0f / %.0f', under_count_1,over_count_1,extinct_1))%[under_count_005],[over_count_005],[extinct_005])
title('1 xMIC')


subplot(2,3,4)
scatter(equiDens_2(:,2),equiDens_2(:,1))
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [5*10^6, 5*10^8];
y = [5*10^6, 5*10^8];
line(x,y,'Color','red','LineWidth',1.5)
text(6*10^6,1.5*10^8,sprintf('%.0f / %.0f / %.0f', under_count_2,over_count_2,extinct_2))%[under_count_005],[over_count_005],[extinct_005])
xlabel('Pf- density (CFU/ml)')
ylabel('Pf+ density (CFU/ml)')
title('2 xMIC')
xlim([5*10^6, 2*10^8])
ylim([5*10^6, 2*10^8])

subplot(2,3,5)
scatter(equiDens_5(:,2),equiDens_5(:,1))
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [1, 10^10];
y = [1, 10^10];
line(x,y,'Color','red','LineWidth',1.5)
text(2,10^9,sprintf('%.0f / %.0f / %.0f', under_count_5,over_count_5,extinct_5))%[under_count_005],[over_count_005],[extinct_005])
xlabel('Pf- density (CFU/ml)')
title('5 xMIC')

subplot(2,3,6)
scatter(equiDens_8(:,2),equiDens_8(:,1))
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [1, 10^10];
y = [1, 10^10];
line(x,y,'Color','red','LineWidth',1.5)
text(2,10^9,sprintf('%.0f / %.0f / %.0f', under_count_8,over_count_8,extinct_8))%[under_count_005],[over_count_005],[extinct_005])
xlabel('Pf- density (CFU/ml)')
title('8 xMIC')


%%

subplot(3,3,7)
scatter(equiDens_15(:,2),equiDens_15(:,1))
hold on;
set(gca,'xscale','log')
set(gca,'yscale','log')
x = [0.5, 5*10^8];
y = [0.5, 5*10^8];
xlim([0.5 10])
ylim([0.5 10])
line(x,y,'Color','red','LineWidth',1.5)




