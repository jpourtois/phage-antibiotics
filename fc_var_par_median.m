function [percent_win,percent_win_n,range2,percent_extinct] = fc_var_par_median(randomScaled, par, phage_start) 

range_dose = 10.^(-2:0.05:2);

if strcmp(par,'theta')
    
    range2 = 0.1:0.025:0.9;
    
elseif strcmp(par,'phi')
    
    range2 = 10.^(5:0.05:7);
    
elseif strcmp(par,'deltaV')
    
    range2 = 10.^(-2:0.125:0);
    
elseif strcmp(par,'epsiR')  
    
    range2 = 0.1:0.025:0.9;
    
elseif strcmp(par,'deltaA') 
    
    range2 = 0.05:0.025:0.5;
    
elseif strcmp(par,'rMin')
    
    range2 = 2:0.5:15;
    
elseif strcmp(par,'init')
    
    range2 = 10.^(2:0.1:8);
    
else 
    
    printf('expression not found')
    
end 


reps = size(randomScaled,1);

n_row = size(range_dose,2);
n_col = size(range2,2);

percent_win = zeros(n_row, n_col);
percent_win_n = zeros(n_row, n_col);
percent_extinct = zeros(n_row, n_col);
percent_extinct_p = zeros(n_row, n_col);
percent_extinct_n = zeros(n_row, n_col);



for k = 1:n_row
    
    dose = range_dose(k);
    
    for j = 1:n_col
        
        sum_win = 0;
        sum_extinct = 0;
        sum_extinct_p = 0;
        sum_extinct_n = 0;
        sum_win_n = 0;
        
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
            %param.aMax = 1; % Peak antibiotic concentration
            param.init = 7*10^7;
            
            %mic = param.kc/(-param.rMin/(param.rMax - param.deltaB) - 1)^(1/param.hill);
            
            
            if strcmp(par,'theta')
                
                param.theta = range2(j);
                
            elseif strcmp(par,'phi')
                
                param.phi = range2(j);
                
            elseif strcmp(par,'deltaV')
                
                param.deltaV = range2(j);
                
            elseif strcmp(par,'epsiR')
                
                param.epsiR = range2(j);
                
            elseif strcmp(par,'deltaA')
                
                param.deltaA = range2(j);
                
            elseif strcmp(par,'rMin')
                
                param.rMin = -range2(j);
                
            elseif strcmp(par,'init')
                
                param.init = range2(j);
                
            end
            
            
            %Set up simulation
            
            total_T = 120; % Total duration of treatment (hour)
            
            nTreatmentPerDay = 2; % Number treatments per day
            %aPerDay = 14*mic; % Total antibiotics dose per day
            
            %param.aMax = aPerDay/nTreatmentPerDay;
            
            param.aMax = dose;
            
            n_treatment = total_T/24*nTreatmentPerDay;
            
            param.MAX_T = total_T/n_treatment;
            totalA = param.aMax*n_treatment;
            
            y0_p = [param.init 10*param.init*phage_start param.aMax];
            y0_n = [param.init 0 param.aMax];
            
            [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param);
            
            y_p_last_day = median(y_p(t_p > 96,1));
            y_n_last_day = median(y_n(t_n > 96,1));
            
            if (y_p_last_day > 1) && (y_p_last_day > y_n_last_day)
                
                sum_win = sum_win + 1;
                
            end
            
            if (y_n_last_day > 1) && (y_n_last_day > y_p_last_day)
                
                sum_win_n = sum_win_n + 1;
                
            end
            
            if (y_p_last_day < 1) && (y_n_last_day < 1)
                
                sum_extinct = sum_extinct + 1;
            end
                
            if (y_p_last_day < 1)
                
                sum_extinct_p = sum_extinct_p + 1;
            end
                
            if (y_n_last_day < 1)  
                
                sum_extinct_n = sum_extinct_n + 1;
                
            end
                
        end
        
        percent_win(k,j) = sum_win/reps;
        
        percent_win_n(k,j) = sum_win_n/reps;
        
        percent_extinct(k,j) = sum_extinct/reps;
        percent_extinct_p(k,j) = sum_extinct_p/reps;
        percent_extinct_n(k,j) = sum_extinct_n/reps;
        
        
        
    end
end




end 