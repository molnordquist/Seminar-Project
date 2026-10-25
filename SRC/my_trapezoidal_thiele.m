%Trapezoidal method
function Vnext = my_trapezoidal_thiele(t,V,h,delta,P,S,mu)
    % One step of the implicit trapezoidal method for Thiele's differential equation
    % Parameters: 
        % t = time 
        % V = reserves at current time 
        % delta = force of interest at time t 
        % P = annual premium rate that is payable at time t
        % S = sum insured payable on death at time t
        % mu = force of mortality

    % Force of Mortality at the current and next time
    mu_curr = mu(t);
    mu_next = mu(t+h);
    % Breaking up Thiele's ODE the Coefficients of V' = a(t)*V + b(t)
    a_n = delta + mu_curr;
    a_next = delta + mu_next;
    b_n = P - mu_curr*S;
    b_next = P - mu_next*S;

    % Implicit trapezoidal 
    numerator = V + (h/2)*(a_n*V + b_n + b_next);
    denominator = 1 - (h/2)*a_next;

    Vnext = numerator/denominator;
end