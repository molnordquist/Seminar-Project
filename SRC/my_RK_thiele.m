%Runge Kutta
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


    %Settig up method
    s1 = f(t,V);
    s2= f(t+h/2, V+(h/2)*s1);
    s3=f(t+h/2, V+(h/2)*s2);
    s4 = f(t+h, V+h*s3);
    Vnext =V + (h/6)*(s1 +2*s2 +2*s3 +s4) ;
 
end