function [sort_pf, sort_no, sort_ratio, p_pf_table, p_no_table, p_ratio_table] = PRCC_function(randomScaled, aMax_rand, aMax)

reps = size(randomScaled,1);
n_par = size(randomScaled,2);

y_p_last_day = zeros(reps,1);
y_n_last_day = zeros(reps,1);
pf_pos_wins = zeros(reps,1);

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
    param.aMax = randomScaled(i,16); % Peak antibiotic concentration
    param.init = randomScaled(i,17); % Starting concentration of both strains
    
    % Turn on or off for fixed antibiotic concentrations
    if (aMax_rand == 0)
   
        param.aMax = aMax;
        
    end
    
    param.init = 10^7;
    
    %Set up simulation
    
    total_T = 120; % Total duration of treatment (hour)
    
    nTreatmentPerDay = 2; % Number treatments per day
    
    n_treatment = total_T/24*nTreatmentPerDay; % Total number of treatments
    
    param.MAX_T = total_T/n_treatment; % Time between between treatments
    %totalA = param.aMax*n_treatment;
    
    y0_p = [param.init 10*param.init param.aMax];
    y0_n = [param.init 0 param.aMax];
    
    [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
    
    y_p_last_day(i) = median(y_p(t_p > 96,1));
    y_n_last_day(i) = median(y_n(t_n > 96,1));
    
    if (y_p_last_day(i) > 1) && (y_n_last_day(i) > 1) && (y_p_last_day(i) > y_n_last_day(i))
        
        pf_pos_wins(i) = 1;
        
    end
    
end

%PCRR with partialcoor

rho_pf = zeros(1,n_par);
rho_no = zeros(1,n_par);
rho_ratio = zeros(1,n_par);

log_y_p = log(y_p_last_day);
log_y_n = log(y_n_last_day);
ratio = log(y_p_last_day./y_n_last_day);

log_param = log(randomScaled);

if (aMax_rand == 0)
   
        log_param(:,16) = log(aMax);
        
end

log_param(:,17) = log(10^7);

for i = 1:n_par
    
    controlled_variables = log_param;
    controlled_variables(:,i) = [];
    
    [rho_pf(i),p_pf(i)] = partialcorr(log_y_p, log_param(:,i), controlled_variables);
    [rho_no(i),p_no(i)] = partialcorr(log_y_n, log_param(:,i), controlled_variables);
    [rho_ratio(i),p_ratio(i)] = partialcorr(ratio, log_param(:,i), controlled_variables);
    
    
end

rho_pf_table = array2table(transpose(rho_pf), 'RowNames', {'H', 'r_{max}', '\Gamma', '\epsilon_R', '\gamma_R', '\epsilon_K', '\xi','\gamma_K','\delta_B','\delta_A', '\lambda', '\delta_V', 'K_d', '\theta', '\phi', 'a_{Max}', 'init'});
rho_no_table = array2table(transpose(rho_no), 'RowNames', {'H', 'r_{max}', '\Gamma', '\epsilon_R', '\gamma_R', '\epsilon_K', '\xi','\gamma_K','\delta_B','\delta_A', '\lambda', '\delta_V', 'K_d', '\theta', '\phi', 'a_{Max}', 'init'});
rho_ratio_table = array2table(transpose(rho_ratio), 'RowNames', {'H', 'r_{max}', '\Gamma', '\epsilon_R', '\gamma_R', '\epsilon_K', '\xi','\gamma_K','\delta_B','\delta_A', '\lambda', '\delta_V', 'K_d', '\theta', '\phi', 'a_{Max}', 'init'});

p_pf_table = array2table(transpose(p_pf), 'RowNames', {'H', 'r_{max}', '\Gamma', '\epsilon_R', '\gamma_R', '\epsilon_K', '\xi','\gamma_K','\delta_B','\delta_A', '\lambda', '\delta_V', 'K_d', '\theta', '\phi', 'a_{Max}', 'init'});
p_no_table = array2table(transpose(p_no), 'RowNames', {'H', 'r_{max}', '\Gamma', '\epsilon_R', '\gamma_R', '\epsilon_K', '\xi','\gamma_K','\delta_B','\delta_A', '\lambda', '\delta_V', 'K_d', '\theta', '\phi', 'a_{Max}', 'init'});
p_ratio_table = array2table(transpose(p_ratio), 'RowNames', {'H', 'r_{max}', '\Gamma', '\epsilon_R', '\gamma_R', '\epsilon_K', '\xi','\gamma_K','\delta_B','\delta_A', '\lambda', '\delta_V', 'K_d', '\theta', '\phi', 'a_{Max}', 'init'});


sort_pf = sortrows(rho_pf_table);
sort_no = sortrows(rho_no_table);
sort_ratio = sortrows(rho_ratio_table);


end