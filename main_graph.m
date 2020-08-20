%%% Main script: "The Impact of Pf Bacteriophages on the Fitness of Pseudomonas aeruginosa: A Mathematical Modeling Approach
% Generate main text figures
% Last modified: August 19th 2020, Julie Pourtois
% 
% First run 'Parameter values' section
% Each figure then has its own section. 

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

%% Figure 4: Antibiotics properties in indirect competition

%%% Figure 4A: Different antibiotic sequestration constants

% Antibiotic sequestration factor
phis = [0 10^5 10^6 10^7];

% Initial conditions for Pf+ strain
B = 5*10^7;
V = 5*10^8;

% Range of antibiotic concentrations
AC = round(10.^(-2:0.25:2.5),2);

param4A = param;

r = zeros(length(AC), length(phis));

for j = 1:length(phis)
    
    for i = 1:length(AC)
        
        % calculate growth rate for each phi and AC
        param4A.phi = phis(j);
        r(i,j) = growthRate(B,V,AC(i),param4A, 1);
    end   
end

% Initial conditions for Pf- strain
V = 0;

param4A = param;

r_n = zeros(length(AC),1);

for i = 1:length(AC)
    
    % calculate growth rate for each AC
    r_n(i) = growthRate(B,V,AC(i),param4A, 0);
end


figure;
subplot(1,3,1)
plot(AC/mic, -r_n)
hold on
plot(AC/mic, -r(:,1))
plot(AC/mic, -r(:,2))
plot(AC/mic, -r(:,3))
plot(AC/mic, -r(:,4))
set(gca,'xscale','log')
xlabel('Antibiotic concentration (xMIC)')
ylabel('Death rate (h^{-1})')
xlim([10^(-1) 10^3])
legend('Pf (-)','\phi = 0', '\phi = 10^{5}', '\phi = 10^{6}', '\phi = 10^{7}','Location','northwest')

%%% Figure 4B: Dependence on growth rate

% Range of dependence on growth rate
epsiR = [0 0.5 1];

% Initial conditions for Pf+ strain
B = 5*10^7;
V = 5*10^8;

% Range of antibiotic concentrations
AC = round(10.^(-2:0.25:2.5),2);

r_p = zeros(length(AC), length(epsiR)); % Pf+
r_n = zeros(length(AC), length(epsiR)); % Pf-

param4B = param;

for j = 1:length(epsiR)
    
    param4B.epsiR = epsiR(j);
    
    for i = 1:length(AC)
        
        % calculate growth rate for each AC
        r_p(i,j) = growthRate(B,V,AC(i),param4B,1); % Pf+ 
        r_n(i,j) = growthRate(B,0,AC(i),param4B,0); % Pf-
    end    
end

subplot(1,3,2)

% Define legend handles
h = zeros(5,1);
h(1) = plot(1,1,'.','Color', [1 0.4 0], 'LineStyle', 'none','MarkerSize', 10 ); hold on;
h(2) = plot(1,1,'.','Color', [0.2 0.6 1], 'LineStyle', 'none','MarkerSize', 10 ); 
h(3) = plot(1,1,'.','Color', [0 0.4 0], 'LineStyle', 'none','MarkerSize', 10 ); 
h(4) = plot(1,1,'Color', 'black', 'LineStyle', '-'); 
h(5) = plot(1,1,'Color', 'black', 'LineStyle', '--'); 

plot(1,1,'.','Color','white','MarkerSize', 12)

plot(AC/mic, -r_n(:,1),'Color', [1 0.4 0], 'LineStyle', '-')
plot(AC/mic, -r_p(:,1),'Color', [1 0.4 0], 'LineStyle', '--')
plot(AC/mic, -r_n(:,2),'Color', [0.2 0.6 1], 'LineStyle', '-')
plot(AC/mic, -r_p(:,2),'Color', [0.2 0.6 1], 'LineStyle', '--')
plot(AC/mic, -r_n(:,3),'Color', [0 0.4 0], 'LineStyle', '-')
plot(AC/mic, -r_p(:,3),'Color', [0 0.4 0], 'LineStyle', '--')
set(gca,'xscale','log')
xlabel('Antibiotic concentration (xMIC)')
ylabel('Death rate (h^{-1})')
xlim([10^(-1) 10^3])
legend(h, '\epsilon_r = 0','\epsilon_r = 0.5','\epsilon_r = 1','Pf (-) ','Pf (+)','Location','northwest')

%%% Figure 4C: Dependence on antibiotic decay rate

decays = [0.05 0.075 0.1 0.15 0.2 0.25]; % Range of decay rate
peak = [0.5*mic 2.5*mic 5*mic]; % Range of peak antibiotic concentration 

total_T = 120; % Total duration of treatment (hour)

nTreatmentPerDay = 1; % Number treatments per day

diff_p_n = zeros(length(decays), length(peak));

param4C = param;

for i = 1:length(decays)
    
    param4C.deltaA = decays(i);
    
    for j = 1:length(peak)
             
        param4C.aMax = peak(j);
        n_treatment = total_T/24*nTreatmentPerDay;
        
        param4C.MAX_T = total_T/n_treatment;
        totalA = param4C.aMax*n_treatment;
        
        y0_p = [7*10^7 7*10^8 param4C.aMax];
        y0_n = [7*10^7 0 param4C.aMax];
        
        [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param4C);
        
        y_p_last_day = median(y_p(t_p > 96,1));
        y_n_last_day = median(y_n(t_n > 96,1));
        
        if (y_p_last_day < 1) && (y_n_last_day < 1)
            diff_p_n(i,j) = 0;
        else
            diff_p_n(i,j) = log10(y_p_last_day/y_n_last_day);
        end
        
    end
end

subplot(1,3,3)
plot(decays, diff_p_n(:,1))
hold on;
plot(decays, diff_p_n(:,2))
plot(decays, diff_p_n(:,3))
xlabel('Rate of antibiotic decay (hour^{-1})')
ylabel('Log10(Pf (+)/Pf (-))')
leg = legend('0.5 xMIC', '2.5 xMIC', '5 xMIC');
title(leg, 'AMC per dose')

%% Figure 5: Growth dynamics of Pf+ and Pf- strains in indirect competition

%%% Figure 5A: Death rates for different metabolic costs

AC = round(10.^(-2:0.25:2.5),2); % Antibiotic concentrations
thetas = [0.2,0.5,0.8]; % Metabolic costs 

B= 5*10^7;
V= 5*10^8;

r_p = zeros(length(AC),length(thetas));
r_n = zeros(length(AC), 1);

param5A = param;

for j = 1:length(thetas)
    
    param5A.theta = thetas(j);
    
    for i = 1:length(AC)
        
        r_p(i,j) = growthRate(B,V,AC(i),param5A,1);
        r_n(i) = growthRate(B,0,AC(i),param5A,0);
        
    end
end

figure;

ha = tight_subplot(2,3,[.15 .1],.1,.1);

%subplot(2,3,1)
axes(ha(1))
plot(AC/mic, -r_n)
hold on
plot(AC/mic, -r_p(:,1))
plot(AC/mic, -r_p(:,2))
plot(AC/mic, -r_p(:,3))
set(gca,'xscale','log')
xlabel('Antibiotic concentration (xMIC)')
ylabel('Death rate (h^{-1})')
legend('Pf (-)', '\theta = 0.2', '\theta = 0.5', '\theta = 0.8', 'Location', 'northwest')
xlim([10^(-1) 10^3])

%%% Figure 5B: Different AMCs per day

total_T = 120; % Total duration of treatment (hour)

nTreatmentPerDay = 1:1:10; % Number treatments per day
aPerDay = [0.5*mic 2.5*mic 5*mic]; % Total antibiotics dose per day % 1

diff_p_n = zeros(length(nTreatmentPerDay), length(aPerDay));

param5B = param;

for i = 1:length(nTreatmentPerDay)
    for j = 1:length(aPerDay)
             
        param5B.aMax = aPerDay(j)/nTreatmentPerDay(i);
        n_treatment = total_T/24*nTreatmentPerDay(i);
        
        param5B.MAX_T = total_T/n_treatment;
        totalA = param5B.aMax*n_treatment;
        
        y0_p = [7*10^7 7*10^8 param5B.aMax];
        y0_n = [7*10^7 0 param5B.aMax];
        
        [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param5B);
        
        y_p_last_day = median(y_p(t_p > 96,1));
        y_n_last_day = median(y_n(t_n > 96,1));
        
        if (y_p_last_day < 1) && (y_n_last_day < 1)
            diff_p_n(i,j) = -10;
        else
            diff_p_n(i,j) = log10(y_p_last_day/y_n_last_day);
        end
        
    end
end

axes(ha(2))
plot(nTreatmentPerDay, diff_p_n(:,1))
hold on
plot(nTreatmentPerDay, diff_p_n(:,2))
plot(nTreatmentPerDay, diff_p_n(:,3))
%set(gca,'yscale','log')
ylabel('Log_{10}  Pf(+)/Pf(-)')
xlabel('# doses (day^{-1})')
leg = legend('0.5 xMIC', '2.5 xMIC', '5 xMIC');
title(leg, 'AMC per day')
xlim([1 10])

%%% Figure 5C: Different AMC per dose

total_T = 120; % Total duration of treatment (hour)

nTreatmentPerDay = 1:10; % Number treatments per day
peak = [0.5*mic 2.5*mic 5*mic]; % Total antibiotics dose per day 

diff_p_n = zeros(length(nTreatmentPerDay), length(peak));

param5C = param;

for i = 1:length(nTreatmentPerDay)
    for j = 1:length(peak)
             
        param5C.aMax = peak(j);
        n_treatment = total_T/24*nTreatmentPerDay(i);
        
        param5C.MAX_T = total_T/n_treatment;
        totalA = param5C.aMax*n_treatment;
        
        y0_p = [7*10^7 7*10^8 param5C.aMax];
        y0_n = [7*10^7 0 param5C.aMax];
        
        [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param5C);
        
        y_p_last_day = median(y_p(t_p > 96,1));
        y_n_last_day = median(y_n(t_n > 96,1));
        
        if (y_p_last_day < 1) && (y_n_last_day < 1)
            diff_p_n(i,j) = 0; 
        else
            diff_p_n(i,j) = log10(y_p_last_day/y_n_last_day);
        end
        
    end
end

axes(ha(3))
plot(nTreatmentPerDay, diff_p_n(:,1))
hold on
plot(nTreatmentPerDay, diff_p_n(:,2))
plot(nTreatmentPerDay, diff_p_n(:,3))
ylabel('Log_{10}  Pf(+)/Pf(-)')
xlabel('# doses (day^{-1})')
leg = legend('0.5 xMIC', '2.5 xMIC', '5 xMIC');
title(leg, 'AMC per dose')
xlim([1 10])

%%% Figure 5D: Three scenarios

% Both strains die

total_T = 24; % Total duration of treatment (hour)

nTreatmentPerDay = 1; % Number treatments per day
aPerDay = 15*mic ; % Total antibiotics dose per day % 1

param5D = param;

param5D.aMax = aPerDay/nTreatmentPerDay;
n_treatment = total_T/24*nTreatmentPerDay;
        
param5D.MAX_T = total_T/n_treatment;
totalA = param5D.aMax*n_treatment;

y0_p = [7*10^7 7*10^8 param5D.aMax];
y0_n = [7*10^7 0 param5D.aMax];
        
[t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param5D);

axes(ha(6))
plot(t_p, y_p(:,1))
hold on
plot(t_n, y_n(:,1))
xlabel('Time (h)')
ylabel('Bacterial density (CFU/ml)')
set(gca,'yscale','log')
legend('Pf (+)','Pf (-)','Location', 'southeast')
title('15 (xMIC/day)')

% Only Pf+ strain survives

total_T = 120; % Total duration of treatment (hour)

nTreatmentPerDay = 1; % Number treatments per day
aPerDay = 7*mic; % Total antibiotics dose per day 

param5D.aMax = aPerDay/nTreatmentPerDay;
n_treatment = total_T/24*nTreatmentPerDay;
        
param5D.MAX_T = total_T/n_treatment;
totalA = param5D.aMax*n_treatment;

y0_p = [7*10^7 7*10^8 param5D.aMax];
y0_n = [7*10^7 0 param5D.aMax];
        
[t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param5D);

axes(ha(5))
plot(t_p, y_p(:,1))
hold on
plot(t_n, y_n(:,1))
xlabel('Time (h)')
ylabel('Bacterial density (CFU/ml)')
set(gca,'yscale','log')
legend('Pf (+)','Pf (-)','Location', 'southeast')
xlim([0 120])
title('7 (xMIC/day)')

% Both strains survive

total_T = 120; % Total duration of treatment (hour)

nTreatmentPerDay = 1; % Number treatments per day
aPerDay = mic ; % Total antibiotics dose per day % 1

param5D.aMax = aPerDay/nTreatmentPerDay;
n_treatment = total_T/24*nTreatmentPerDay;
        
param5D.MAX_T = total_T/n_treatment;
totalA = param5D.aMax*n_treatment;

y0_p = [7*10^7 7*10^8 param5D.aMax];
y0_n = [7*10^7 0 param5D.aMax];
        
[t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param5D);

axes(ha(4))
plot(t_p, y_p(:,1))
hold on
plot(t_n, y_n(:,1))
set(gca,'yscale','log')
ylim([10^6 10^8])
xlabel('Time (h)')
ylabel('Bacterial density (CFU/ml)')
legend('Pf (+)','Pf (-)','Location', 'southeast')
xlim([0 120])
title('1 (xMIC/day)')

%set(gcf, 'Position',  [50, 50, 1000, 500])


%% Figure 3: Growth dynamics of Pf+ and Pf- in direct competition

%%% Figure 3A: Death rate versus AMC - No dependence on growth rate 

AC = round(10.^(-2:0.5:2.5),2);
B= 5*10^6;
V = 5*10^7;

r_p = zeros(length(AC),1); % Pf+
r_n = zeros(length(AC),1); % Pf-

param3A = param;
param3A.epsiR = 0;

for i = 1:length(AC)
    
    r_p(i) = growthRate(B,V,AC(i),param3A,1);
    r_n(i) = growthRate(B,V,AC(i),param3A,0);
end

figure;
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

%%% Figure 3D: Three different scenarios

% Both strains survive

total_T = 100; % Total duration of treatment (hour)
nTreatmentPerDay = 1; % Number treatments per day
aPerDay = mic; % Total antibiotics dose per day 

param3D = param;

param3D.aMax = aPerDay/nTreatmentPerDay;
n_treatment = total_T/24*nTreatmentPerDay;

param3D.MAX_T = total_T/n_treatment; 
totalA = param3D.aMax*n_treatment;

y0 = [10^7 10^7 10^8 param3D.aMax];

[t,y] = competeRegimens(n_treatment, y0, param3D);

subplot(2,3,4)
plot(t,y(:,1))
hold on
plot(t,y(:,2))
set(gca,'yscale','log')
legend('Pf (+)','Pf (-)','Location', 'southwest')
title('1 (xMIC)') % 0.25
ylabel('Bacterial density (CFU/ml)')
ylim([10^3 10^8])
%text(5,2*10^5,'Peak AMC = 0.5 \mug/ml')

% Only Pf+ survive

total_T = 100; % Total duration of treatment (hour)
nTreatmentPerDay = 1; % Number treatments per day
aPerDay = 5*mic; % Total antibiotics dose per day % 1

param3D.aMax = aPerDay/nTreatmentPerDay;
n_treatment = total_T/24*nTreatmentPerDay;

param3D.MAX_T = total_T/n_treatment; 
totalA = param3D.aMax*n_treatment;

y0 = [10^7 10^7 10^8 param3D.aMax];

[t,y] = competeRegimens(n_treatment, y0,param3D);

subplot(2,3,5)
plot(t,y(:,1))
hold on
plot(t,y(:,2))
set(gca,'yscale','log')
legend('Pf (+)','Pf (-)','Location', 'southwest')
xlabel('Time (h)')
title('5 (xMIC)')
ylim([10^3 10^8])

% Both strains die

total_T = 100; % Total duration of treatment (hour)
nTreatmentPerDay = 1; % Number treatments per day
aPerDay = 15*mic; % Total antibiotics dose per day % 1

param3D.aMax = aPerDay/nTreatmentPerDay;
n_treatment = total_T/24*nTreatmentPerDay;

param3D.MAX_T = total_T/n_treatment; 
totalA = param3D.aMax*n_treatment;

y0 = [10^7 10^7 10^8 param3D.aMax];

[t,y] = competeRegimens(n_treatment, y0,param3D);

subplot(2,3,6)
plot(t,y(:,1))
hold on
plot(t,y(:,2))
set(gca,'yscale','log')
legend('Pf (+)','Pf (-)')
title('15 (xMIC)')

%% Figure 2: Growth of pf+ and pf- without antibiotics

%%% Figure 2A: Bacterial density at equilibrium
range_theta = 0:0.05:0.8;

b_equi_mat = zeros(length(range_theta), 1);

K = param.K;
rMax = param.rMax;
deltaB = param.deltaB;

for i = 1:length(range_theta)

    theta = range_theta(i);
 
    b_equi_mat(i) = K*(1-(deltaB/(rMax*(1 - theta)))^(1/param.exponent));

end

%%% Figure 2B: Competition without antibiotics

MAX_T = 3000;
y0 = [5*10^7 100 5*10^8 0];
param2B = param;

[t_10,y_10] = ode45(@(t,y) EmaxcompeteNoResist(t,y,param2B), [0 MAX_T], y0);

param2B.theta = 0.4;
[t_04,y_04] = ode45(@(t,y) EmaxcompeteNoResist(t,y,param2B), [0 MAX_T], y0);

param2B.theta = 0.6;
[t_06,y_06] = ode45(@(t,y) EmaxcompeteNoResist(t,y,param2B), [0 MAX_T], y0);


figure;
subplot(1,2,1)

plot(range_theta, b_equi_mat,'LineWidth',1)
hold on;
xlabel('Metabolic cost','FontSize',12)
ylabel('Equilibrium density (CFU/ml)','FontSize',12)
hold off

subplot(1,2,2)

h = zeros(5,1);
h(1) = plot(t_10/24, y_10(:,1),'Color',[0, 0.4470, 0.7410],'LineWidth',1); hold on;
h(2) = plot(t_04/24, y_04(:,1),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1);
h(3) = plot(t_06/24, y_06(:,1),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',1);
h(4) = plot(t_06/24, y_06(:,1),'Color','black','LineWidth',1);
h(5) = plot(t_06/24, y_06(:,2),'Color','black','LineWidth',1,'LineStyle','--');

plot(t_10/24, y_10(:,1),'Color',[0, 0.4470, 0.7410],'LineWidth',1)
plot(t_10/24, y_10(:,2),'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',1)
plot(t_04/24, y_04(:,1),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1)
plot(t_04/24, y_04(:,2),'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',1)
plot(t_06/24, y_06(:,1),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',1)
plot(t_06/24, y_06(:,2),'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--','LineWidth',1)
set(gca,'yscale','log')
ylabel('Bacterial density (CFU/ml)','FontSize',12)
xlabel('Time (days)','FontSize',12)
legend(h, '\theta = 0.2', '\theta = 0.4','\theta = 0.6','Pf (+)','Pf (-)','location','southeast')
xlim([0 120])
ylim([1 10^8])

set(gcf, 'Position',  [50, 50, 800, 300])



