function [t_p, y_p, t_n, y_n] = compareRegimens(n_treatment, y0_p, y0_n, param)

t_p = 0;
t_n = 0;
y_p = zeros(1, 3);
y_n = zeros(1, 3);

for k = 1:n_treatment
    
    % Phage-positive
    if (k == 1)
        y0 = y0_p;
    else
        y0 = y_p(length(y_p(:,1)),:);
    end
    
    y0(3) = param.aMax;
    
    param.pf = 1;
    
    [t,y] = ode45(@(t,y) EmaxcompareNoResist(t,y,param), [0 param.MAX_T], y0) ;
    
    t_p = [t_p; (t+t_p(end))];
    y_p = [y_p; y];
    
    % Phage-negative
    
    if (k == 1)
        y0 = y0_n;
    else
        y0 = y_n(length(y_n(:,1)),:);
    end
    
    y0(3) = param.aMax;
    
    param.pf = 0;
    
    [t,y] = ode45(@(t,y) EmaxcompareNoResist(t,y,param), [0 param.MAX_T], y0);
    
    t_n = [t_n; (t+t_n(end))];
    y_n = [y_n; y];
    
end 