function dV = thieles(t,V,delta,P,S,mu)
    % Thiele's ODE modelnin terms of constant mortality model 
    % Parameters: 
        % dv = change in reserves over time 
        % t = time 
        % V = reserves at current time 
        % delta = force of interest at time t 
        % P = annual premium rate that is payable at time t
        % S = sum insured payable on death at time t
        % mu = force of mortality
    dV = delta*V + P - mu*(S-V);
end