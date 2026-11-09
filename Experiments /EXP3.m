
%% Exact Analytical Solution for Gompertz Mortality

% Parameters from previous textbook example
T = 20;
S = 100000;
P = 2500;
delta = 0.04;
x = 30;

% Seeting up Gompertz mortality parameters
B = 2.7e-6;
c = 1.124;
mu = @(t) B*c.^(x+t);

% Creating derived Integrating factor
I = @(t) exp(-delta*t - (B*c.^(x+t))/log(c));

% Creating time step grid
t = linspace(0,T,401);

% Initializing reserve array
V_exact_gompertz = zeros(1,length(t));

% Creating function for integral
int_component = @(s) I(s).*(P-S*mu(s));
% Calculate analytical reserve at each time
for n = 1:length(t)
    tn = t(n);
    % Calculating the integral from xmin=tn to xmax=T
    integral_sol = integral(int_component,tn,T);
    V_exact_gompertz(n) = (I(T)*S - integral_sol)/I(tn);

end

% Find reserve at year 10, indicated year of textbook. Same as 40 years old
YR10 = find(t == 10,1);

fprintf('Gompertz reserve at year 10: %.10f\n',V_exact_gompertz(YR10));

