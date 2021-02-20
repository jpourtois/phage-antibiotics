% Figure 2 A,B
% Julie Pourtois
% Last modified: Feb 19th 2021

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
