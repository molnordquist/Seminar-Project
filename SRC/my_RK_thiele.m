%Runge-Kutta
function Vnext = my_RK_thiele(t,V,h,delta,P,S,mu)
    % One step of RK4's method for Thiele's differential equation
    % Parameters: 
        % t = time 
        % V = reserves at current time 
        % delta = force of interest at time t 
        % P = annual premium rate that is payable at time t
        % S = sum insured payable on death at time t
        % mu = force of mortality
        
    f = @(t, V) delta*V + P - mu(t)*(S-V);


    %Settig up four slope estimate
    k1 = f(t,V);
    k2= f(t+h/2, V+(h/2)*k1);
    k3=f(t+h/2, V+(h/2)*k2);
    k4 = f(t+h, V+h*k3);
    Vnext =V + (h/6)*(k1 +2*k2 +2*k3 +k4) ; %weighted average
 
end