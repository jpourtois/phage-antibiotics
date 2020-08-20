%%% Sensitivity analysis - Theta,phi,Kd
% Generate supplemental figures S3
% Last modified: August 19th 2020, Julie Pourtois

%%
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
param.phi = 10^6; % Antibiotic sequestration factor
param.aMax = 1; % Peak antibiotic concentration


pMin(1) = 0.1; % Theta 
pMax(1) = 0.9;

pMin(2) = 10^5; % Phi
pMax(2) = 10^7;

pMin(3) = 10^(13); % Kd
pMax(3) = 10^(15);

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
randomUnscaled = lhsdesign(reps,3);

% Scale the random number with the range spread in log space
randomdiffLog = zeros(reps,3);
 
for n = 1:3
    randomdiffLog(:,n) = randomUnscaled(:,n)*diffLog(n);    
end

randomLogScaled = logMin + randomdiffLog;

%Get out of log space
randomScaled = exp(randomLogScaled);

equiDens_005 = zeros(reps,6);
equiDens_05 = zeros(reps,6);
equiDens_1 = zeros(reps,6);
equiDens_2 = zeros(reps,6);
equiDens_5 = zeros(reps,6);
equiDens_8 = zeros(reps,6);
equiDens_15 = zeros(reps,6);

over_count_005 = 1;
over_count_05 = 1;
over_count_1 = 1;
over_count_2 = 1;
over_count_5 = 1;
over_count_8 = 1;
over_count_15 = 1;

for i = 1:reps
    
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
    
    % Change theta
    param.theta = randomScaled(i,1);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_005(i,1) = y_p_last_day;
    equiDens_005(i,2) = y_n_last_day;
    
    param.theta = 0.2;
    
    % Change phi
    
    param.phi = randomScaled(i,2);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_005(i,3) = y_p_last_day;
    equiDens_005(i,4) = y_n_last_day;
    
    param.phi = 10^6;
    
    % Change Kd
    
    param.kD = randomScaled(i,3);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_005(i,5) = y_p_last_day;
    equiDens_005(i,6) = y_n_last_day;
    
    param.kD = 10^(14);

    % 0.5 xMIC
   
    param.aMax = 0.5*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    % Change theta
    param.theta = randomScaled(i,1);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_05(i,1) = y_p_last_day;
    equiDens_05(i,2) = y_n_last_day;
    
    param.theta = 0.2;
    
    % Change phi
    
    param.phi = randomScaled(i,2);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_05(i,3) = y_p_last_day;
    equiDens_05(i,4) = y_n_last_day;
    
    param.phi = 10^6;
    
    % Change Kd
    
    param.kD = randomScaled(i,3);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_05(i,5) = y_p_last_day;
    equiDens_05(i,6) = y_n_last_day;
    
    param.kD = 10^(14);
    
    % 1 xMIC
    
    param.aMax = mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    % Change theta
    param.theta = randomScaled(i,1);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_1(i,1) = y_p_last_day;
    equiDens_1(i,2) = y_n_last_day;
    
    param.theta = 0.2;
    
    % Change phi
    
    param.phi = randomScaled(i,2);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_1(i,3) = y_p_last_day;
    equiDens_1(i,4) = y_n_last_day;
    
    param.phi = 10^6;
    
    % Change Kd
    
    param.kD = randomScaled(i,3);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_1(i,5) = y_p_last_day;
    equiDens_1(i,6) = y_n_last_day;
    
    param.kD = 10^(14);
    
    % 2 xMIC
  
    param.aMax = 2*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    % Change theta
    param.theta = randomScaled(i,1);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_2(i,1) = y_p_last_day;
    equiDens_2(i,2) = y_n_last_day;
    
    param.theta = 0.2;
    
    % Change phi
    
    param.phi = randomScaled(i,2);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_2(i,3) = y_p_last_day;
    equiDens_2(i,4) = y_n_last_day;
    
    param.phi = 10^6;
    
    % Change Kd
    
    param.kD = randomScaled(i,3);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_2(i,5) = y_p_last_day;
    equiDens_2(i,6) = y_n_last_day;
    
    param.kD = 10^(14);
    
    % 5 xMIC
    
    param.aMax = 5*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    % Change theta
    param.theta = randomScaled(i,1);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_5(i,1) = y_p_last_day;
    equiDens_5(i,2) = y_n_last_day;
    
    param.theta = 0.2;
    
    % Change phi
    
    param.phi = randomScaled(i,2);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_5(i,3) = y_p_last_day;
    equiDens_5(i,4) = y_n_last_day;
    
    param.phi = 10^6;
    
    % Change Kd
    
    param.kD = randomScaled(i,3);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_5(i,5) = y_p_last_day;
    equiDens_5(i,6) = y_n_last_day;
    
    param.kD = 10^(14);
    
    % 8 xMIC
    param.aMax = 8*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    % Change theta
    param.theta = randomScaled(i,1);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_8(i,1) = y_p_last_day;
    equiDens_8(i,2) = y_n_last_day;
    
    param.theta = 0.2;
    
    % Change phi
    
    param.phi = randomScaled(i,2);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_8(i,3) = y_p_last_day;
    equiDens_8(i,4) = y_n_last_day;
    
    param.phi = 10^6;
    
    % Change Kd
    
    param.kD = randomScaled(i,3);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_8(i,5) = y_p_last_day;
    equiDens_8(i,6) = y_n_last_day;
    
    param.kD = 10^(14);
    
    % 15 xMIC
    
    param.aMax = 15*mic;
    n_treatment = total_T/24*nTreatmentPerDay;
    
    param.MAX_T = total_T/n_treatment;
    
    y0_p = [7*10^7 7*10^8 param.aMax];
    y0_n = [7*10^7 0 param.aMax];
    
    % Change theta
    param.theta = randomScaled(i,1);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_15(i,1) = y_p_last_day;
    equiDens_15(i,2) = y_n_last_day;
    
    param.theta = 0.2;
    
    % Change phi
    
    param.phi = randomScaled(i,2);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_15(i,3) = y_p_last_day;
    equiDens_15(i,4) = y_n_last_day;
    
    param.phi = 10^6;
    
    % Change Kd
    
    param.kD = randomScaled(i,3);
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day = median(y_p(t_p > 96,1));
    y_n_last_day = median(y_n(t_n > 96,1));
    
    no_extinct = y_p_last_day > 1 & y_n_last_day > 1;
    
    equiDens_15(i,5) = y_p_last_day;
    equiDens_15(i,6) = y_n_last_day;
    
    param.kD = 10^(14);
     
end


%%
figure()
data = {equiDens_005(:,[1 3 5 6]),equiDens_05(:,[1 3 5 6]),equiDens_1(:,[1 3 5 6]),equiDens_2(:,[1 3 5 6]), equiDens_5(:,[1 3 5 6]),equiDens_8(:,[1 3 5 6]),equiDens_15(:,[1 3 5 6])};
boxplotGroup(data,'PrimaryLabels', {'0.05', '0.5','1','2','5','8','15'},'SecondaryLabels',{'\Theta (0.1-0.9)' '\Phi (10^{5}-10^{7})' 'K_{d} (10^{13}-10^{15})' 'Pf-'})
xlabel({'','Antibiotic concentration (xMIC)'},'fontsize',12)
ylabel('Bacterial density (CFU/ml)','fontsize',12)

set(gcf, 'Position',  [200, 200, 800, 500])
