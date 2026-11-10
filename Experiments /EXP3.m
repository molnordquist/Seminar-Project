clear;
clc;
close all;
addpath('../SRC');

%% Analytical Solution for Constant Mortality
% Parameters from previous textbook example
T = 20;
S = 100000;
P = 2500;
delta = 0.04;
x = 30; %age

% Setting up Gompertz mortality parameters
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

fprintf('Constant reserve at year 10: %.10f\n',V_constant(YR10));

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
YR10 = find(t == 10,1);

fprintf('Gompertz reserve at year 10: %.10f\n',V_exact_gompertz(YR10));

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
YR10 = find(t == 10,1);

fprintf('Makeham reserve at year 10: %.10f\n',V_exact_makeham(YR10));

%% Plot all three analytical solutions
figure;
plot(t, V_constant, 'LineWidth',2);
hold on;
plot(t, V_exact_makeham, 'LineWidth',2);
plot(t, V_exact_gompertz, 'LineWidth',2);
xlabel('Time (yrs)');
ylabel('Reserve Amount');
title('Analytical Solutions for Mortality Laws in Thieles ODE');
legend('Constant Mortality','Gompertz','Makeham','Location','best');
grid on;
hold off;

%% Convergence Analysis 

%Define parameters
T = 20; 
delta = 0.04; 
P = 2500; 
S = 100000; 

mu_const = 0.01; 
mu = @(t) mu_const; %mu is set up as a function of t


% Test h = -0.05 first (solving backwards)
N=400;
h= -T/N;
%creating equally spaced time-steps
t = linspace(T,0,N+1); 

%Initializing arrays
V_euler =zeros(1, N+1);
V_RK = zeros(1, N+1);
V_trap = zeros(1, N+1);

% Setting the first term in the array as the Terminal condition 
V_euler(1) = 100000;
V_trap(1) = 100000;
V_RK(1) = 100000;

%running methods 
for n = 1:N
    V_euler(n+1) = my_Euler_thiele(t(n),V_euler(n),h,delta,P,S,mu);
    V_trap(n+1) = my_trapezoidal_thiele(t(n),V_trap(n),h,delta,P,S,mu);
    V_RK(n+1) = my_RK_thiele(t(n),V_RK(n),h,delta,P,S,mu);
end

% Compare reserves at year 10
% finding spot in array of year 10
YR10 = find(t==10,1);
%using coded constant exact solution
V_exact_10 = V_exact_constant(10);
error_euler = abs(V_euler(YR10)-V_exact_10);
error_RK = abs(V_RK(YR10)-V_exact_10);
error_trap = abs(V_trap(YR10)-V_exact_10);

%Reserve values by method
fprintf('Step size = %.4f\n', h)
fprintf('Error of Euler at year 10: %.10f\n', error_euler);
fprintf('Error of RK4 at year 10: %.10f\n', error_RK);
fprintf('Error of Trapezoidal at year 10: %.10f\n', error_trap);

%%
% New step-size -0.025

N=800;
h2= -T/N;
%creating equally spaced time-steps
t = linspace(T,0,N+1); 

%Initializing arrays
V_euler2 =zeros(1, N+1);
V_RK2 = zeros(1, N+1);
V_trap2 = zeros(1, N+1);

% Setting the first term in the array as the Terminal condition 
V_euler2(1) = 100000;
V_trap2(1) = 100000;
V_RK2(1) = 100000;

%running methods 
for n = 1:N
    V_euler2(n+1) = my_Euler_thiele(t(n),V_euler2(n),h2,delta,P,S,mu);
    V_trap2(n+1) = my_trapezoidal_thiele(t(n),V_trap2(n),h2,delta,P,S,mu);
    V_RK2(n+1) = my_RK_thiele(t(n),V_RK2(n),h2,delta,P,S,mu);
end

% Compare reserves at year 10
% finding spot in array of year 10
YR10 = find(t==10,1);
%using coded constant exact solution
V_exact_10 = V_exact_constant(10);
error_euler2 = abs(V_euler2(YR10)-V_exact_10);
error_RK2 = abs(V_RK2(YR10)-V_exact_10);
error_trap2 = abs(V_trap2(YR10)-V_exact_10);

%Reserve values by method
fprintf('Step size = %.4f\n', h2)
fprintf('Error of Euler at year 10: %.10f\n', error_euler2);
fprintf('Error of RK4 at year 10: %.10f\n', error_RK2);
fprintf('Error of Trapezoidal at year 10: %.10f\n', error_trap2);
%%
% New step-size -0.025

N=1600;
h3= -T/N;
%creating equally spaced time-steps
t = linspace(T,0,N+1); 

%Initializing arrays
V_euler3 =zeros(1, N+1);
V_RK3 = zeros(1, N+1);
V_trap3 = zeros(1, N+1);

% Setting the first term in the array as the Terminal condition 
V_euler3(1) = 100000;
V_trap3(1) = 100000;
V_RK3(1) = 100000;

%running methods 
for n = 1:N
    V_euler3(n+1) = my_Euler_thiele(t(n),V_euler3(n),h3,delta,P,S,mu);
    V_trap3(n+1) = my_trapezoidal_thiele(t(n),V_trap3(n),h3,delta,P,S,mu);
    V_RK3(n+1) = my_RK_thiele(t(n),V_RK3(n),h3,delta,P,S,mu);
end

% Compare reserves at year 10
% finding spot in array of year 10
YR10 = find(t==10,1);
%using coded constant exact solution
V_exact_10 = V_exact_constant(10);
error_euler3 = abs(V_euler3(YR10)-V_exact_10);
error_RK3 = abs(V_RK3(YR10)-V_exact_10);
error_trap3 = abs(V_trap3(YR10)-V_exact_10);

%Reserve values by method
fprintf('Step size = %.4f\n', h3)
fprintf('Error of Euler at year 10: %.10f\n', error_euler3);
fprintf('Error of RK4 at year 10: %.10f\n', error_RK3);
fprintf('Error of Trapezoidal at year 10: %.10f\n', error_trap3);

%%
%calculating order of convergence 
k_euler = abs(log(error_euler3/error_euler2)/log(2));
fprintf('Euler order of convergence: %.4f\n',k_euler);
k_trap = abs(log(error_trap3/error_trap2)/log(2));
fprintf('Trapezoidal order of convergence: %.4f\n',k_trap);
k_RK = abs(log(error_RK3/error_RK2)/log(2));
fprintf('RK4 order of convergence: %.4f\n',k_RK);

