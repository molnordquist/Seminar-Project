clear;
clc;
close all;
addpath('../SRC');

exportgraphics(gcf, fullfile(figureFolder, 'Convergence_Plot.png'),'Resolution', 300);

%% Analytical Solution for Constant Mortality
% Parameters from previous textbook example
T = 20;
S = 100000;
P = 2500;
delta = 0.04;
x = 30; %age

% mortality
mu_const = 0.01;

%creating timegrid
t=linspace(0,T,401);

%Setting up components
A= delta +mu_const;
B= P-mu_const*S;

% Creating exact reserve
V_exact_constant = @(t) -B/A + (S + B/A)*exp(A*(t-T));
V_constant = V_exact_constant(t);

% Find reserve at year 10, indicated year of textbook. Same as 40 years old
YR10 = find(t == 10,1);
YR20 = find(t == T,1);
YR0 = find(t == 0,1);
fprintf('Constant reserve at year 20: %.10f\n',V_constant(YR20));
fprintf('Constant reserve at year 10: %.10f\n',V_constant(YR10));
fprintf('Constant reserve at year 0: %.10f\n',V_constant(YR0));
%% Analytical Solution for Gompertz Mortality

% Parameters from previous textbook example
T = 20;
S = 100000;
P = 2500;
delta = 0.04;
x = 30;

% Setting up Gompertz mortality parameters
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
full_int_gomp = @(s) I(s).*(P-S*mu(s));
% Calculate analytical reserve at each time
for n = 1:length(t)
    tn = t(n);
    % Calculating the integral from xmin=tn to xmax=T
    integral_sol_gomp = integral(full_int_gomp,tn,T);
    V_exact_gompertz(n) = (I(T)*S - integral_sol_gomp)/I(tn);
end

% Find reserve at year 10, indicated year of textbook. Same as 40 years old
YR0 = find(t == 0,1);
YR10 = find(t == 10,1);
YR20 = find(t == T,1);
fprintf('Gompertz reserve at year 0: %.10f\n',V_exact_gompertz(YR0));
fprintf('Gompertz reserve at year 10: %.10f\n',V_exact_gompertz(YR10));
fprintf('Gompertz reserve at year 20: %.10f\n',V_exact_gompertz(YR20));


%%
%Analytical Solution for Makeham's
% Parameters from previous textbook example
T = 20;
S = 100000;
P = 2500;
delta = 0.04;
x = 30;

% Seeting up Gompertz mortality parameters
A = 0.00022;
B = 2.7e-6;
c = 1.124;
mu = @(t) A + B*c.^(x+t);

% Creating derived Integrating factor
I = @(t) exp(-delta*t -A*t - (B*c.^(x+t))/log(c));

% Creating time step grid
t = linspace(0,T,401);

% Initializing reserve array
V_exact_makeham = zeros(1,length(t));


% Creating function for integral
full_int_mak = @(s) I(s).*(P-S*mu(s));
% Calculate analytical reserve at each time
for n = 1:length(t)
    tn = t(n);
    % Calculating the integral from xmin=tn to xmax=T
    integral_sol_mak = integral(full_int_mak,tn,T);
    V_exact_makeham(n) = (I(T)*S - integral_sol_mak)/I(tn);
end

% Find reserve at year 10, indicated year of textbook. Same as 40 years old
YR0 = find(t == 0,1);
YR10 = find(t == 10,1);
YR20 = find(t == T,1);

fprintf('Makeham reserve at year 0: %.10f\n',V_exact_makeham(YR0));
fprintf('Makeham reserve at year 10: %.10f\n',V_exact_makeham(YR10));
fprintf('Makeham reserve at year 20: %.10f\n',V_exact_makeham(YR20));

%%
%check using terminal condition 
fprintf('Constant mortality error: %.10f\n',abs(V_constant(YR20) -S));
fprintf('Gompertz mortality error: %.10f\n',abs(V_exact_gompertz(YR20) -S));
fprintf('Makeham mortality error: %.10f\n',abs(V_exact_makeham(YR20) -S));
%% Plot all three analytical solutions
figure;
plot(t, V_constant, 'LineWidth',2);
hold on;
plot(t, V_exact_makeham, 'LineWidth',2);
plot(t, V_exact_gompertz, 'LineWidth',2);
xlabel('Time (yrs)');
ylabel('Reserve Amount');
title('Analytical Solutions for Mortality Laws in Thieles ODE');
legend('Constant Mortality','Makeham','Gompertz','Location','best');
grid on;
hold off;

