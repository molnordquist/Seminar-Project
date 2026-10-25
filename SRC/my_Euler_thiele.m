%Euler's method 
function Vnext = my_Euler_thiele(V,h,delta,P,S,mu_curr)
    % One step of Euler's method for Thiele's differential equation
    % Parameters: 
        % t = time 
        % V = reserves at current time 
        % delta = force of interest at time t 
        % P = annual premium rate that is payable at time t
        % S = sum insured payable on death at time t
        % mu = force of mortality

    % Force of Mortality at the current and next time
    mu_curr = mu(t);

    %Thieles
    dV = delta*V + P - mu_curr*(S-V);

    % Breaking up Thiele's ODE the Coefficients of V' = a(t)*V + b(t)
    Vnext = V + h*dV;
 
end