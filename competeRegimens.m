function [t_p, y_p] = competeRegimens(n_treatment,y_init,param)

t_p = 0;
y_p = zeros(1, 4);

for k = 1:n_treatment
    
    
    if (k == 1)
        y0 = y_init;
    else
        y0 = y_p(length(y_p(:,1)),:);
    end
    
    y0(4) = param.aMax;
    
    [t,y] = ode45(@(t,y) EmaxcompeteNoResist(t,y,param), [0 param.MAX_T], y0) ;
    
    t_p = [t_p; (t+t_p(end))];
    y_p = [y_p; y];
    
    
end 
