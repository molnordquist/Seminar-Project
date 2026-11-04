
% Textbook ex 7.13 from Actuarial Mathematics for Life Contingent Risks
clear;
clc;
close all;
addpath('../SRC');

%%
%Error using exact analytical solution
T = 20;
S = 100000;
P = 2500;
delta = 0.04;
mu_const = 0.01;
mu = @(t) mu_const; %randomly chosen but can be estimated through mortality data 
%setting up equation
a = delta + mu_const;
b = P - mu_const*S;
V_exact = @(t) -b/a + (S+b/a)*exp(a*(t-T));
% Time grid
N = 400;
h = -T/N;
t = linspace(T,0,N+1);

% Initialize
V_euler = zeros(1,N+1);
V_trap = zeros(1,N+1);
V_RK4 = zeros(1,N+1);
%Setting first entry in array as terminal condition
V_euler(1) = S;
V_trap(1) = S;
V_RK4(1) = S;

% Numerical approximations
tic
for n = 1:N
    V_euler(n+1) = my_Euler_thiele(t(n),V_euler(n),h,delta,P,S,mu);
end
compTime_euler05 = toc;

tic
for n = 1:N
V_trap(n+1) = my_trapezoidal_thiele(t(n),V_trap(n),h,delta,P,S,mu);
end 
compTime_Trap05 = toc;

tic
for n = 1:N
    V_RK4(n+1) = my_RK_thiele( t(n),V_RK4(n),h,delta,P,S,mu);
end
compTime_RK05 = toc;
%Print computational time results
fprintf('\nComputational Time Comparison\n');
fprintf('Euler: %.8f seconds\n',compTime_euler05);
fprintf('Trapezoidal: %.8f seconds\n',compTime_Trap05);
fprintf('RK4: %.8f seconds\n',compTime_RK05);

% Find year 10
YR10 = find(t == 10,1);

% Exact reserve
V_true = V_exact(10);

% True numerical errors
error_euler = abs(V_euler(YR10)-V_true);
error_trap = abs(V_trap(YR10)-V_true);
error_RK = abs(V_RK4(YR10)-V_true);

fprintf('Exact reserve: %.10f\n',V_true);
fprintf('Euler error: %.10f\n',error_euler);
fprintf('Trapezoidal error: %.10f\n',error_trap);
fprintf('RK4 error: %.10f\n',error_RK);
%% 
% smaller step size
% Makeham mortality with smaller step size h = -0.01

T = 20;
S = 100000;
P = 2500;
delta = 0.04;
mu = @(t) 0.01;

% New step size
N2 = 2000;
h2 = -T/N2;

% New time grid
t2 = linspace(T,0,N2+1);

% Initializing the arrays
V_euler01= zeros(1,N2+1);
V_trap01 = zeros(1,N2+1);
V_RK01 = zeros(1,N2+1);

% Terminal condition
V_euler01(1) = S;
V_trap01(1) = S;
V_RK01(1) = S;

% Run numerical methods
tic
for n = 1:N2
    V_euler01(n+1) = my_Euler_thiele(t2(n),V_euler01(n),h2,delta,P,S,mu);
end
compTime_euler01 = toc;

tic 
for n=1:N2
    V_trap01(n+1) = my_trapezoidal_thiele(t2(n),V_trap01(n),h2,delta,P,S,mu);
end
compTime_Trap01 = toc;

tic
for n=1:N2
    V_RK01(n+1) = my_RK_thiele(t2(n),V_RK01(n),h2,delta,P,S,mu);
end
compTime_RK01 = toc;
%Print computational time results
fprintf('\nComputational Time Comparison\n');
fprintf('Euler: %.8f seconds\n',compTime_euler01);
fprintf('Trapezoidal: %.8f seconds\n',compTime_Trap01);
fprintf('RK4: %.8f seconds\n',compTime_RK01);

% Find year 10
YR10_2 = find(t2 == 10, 1);
% Exact reserve
V_true = V_exact(10);

% True numerical errors
error_euler01 = abs(V_euler01(YR10_2)-V_true);
error_trap01 = abs(V_trap01(YR10_2)-V_true);
error_RK01 = abs(V_RK01(YR10_2)-V_true);

fprintf('Exact reserve: %.10f\n',V_true);
fprintf('Euler error: %.10f\n',error_euler01);
fprintf('Trapezoidal error: %.10f\n',error_trap01);
fprintf('RK4 error: %.10f\n',error_RK01);


%%
%comparing errors for both step sizes 
%difference in computational time for step-sizes
diffComp_eul = abs(compTime_euler05 - compTime_euler01);
diffComp_trap= abs(compTime_Trap05 - compTime_Trap01);
diffComp_RK = abs(compTime_RK05 - compTime_RK01);
%print 
fprintf('\nDifference in Computational Time (h=-0.05 vs h=-0.01)\n');
fprintf('Euler: %.8f seconds\n',diffComp_eul);
fprintf('Trapezoidal: %.8f seconds\n',diffComp_trap);
fprintf('RK4: %.8f seconds\n',diffComp_RK);

% difference in errors for both step sizes
diffError_eul = abs(error_euler -error_euler01);
diffError_trap = abs(error_trap -error_trap01);
diffError_RK = abs(error_RK -error_RK01);
%print
fprintf('\nDifference in Numerical Errors (h=-0.05 vs h=-0.01)\n');
fprintf('Euler: %.10f\n',diffError_eul);
fprintf('Trapezoidal: %.10f\n',diffError_trap);
fprintf('RK4: %.10f\n',diffError_RK);

%%
% Additional analysis
%% 
%Define parameters
T = 20; %policy term is 20 years
delta = 0.04; %interest rate per year
P = 2500; %Policy holder pays continuous premiums
S = 100000; %death benefit
% we are using Standard Makeham mortality for this example
mu = @(t) .00022 +(2.7e-6)*1.124.^(30+t); %constant component and age-dependent component

% Test h = -0.05 first (solving backwards)
N=400;
h= -T/N;
%creating equally spaced time-steps
t = linspace(T,0,N+1); %start condition T, end condition 0, N+1 points
% t=[20, 19.95, 19.90, ... , 0]

% Setting the first term in the array as the Terminal condition = endowment benefit
V_euler_M(1) = 100000;
V_trap_M(1) = 100000;
V_RK_M(1) = 100000;

%running methods 
for n = 1:N
    V_euler_M(n+1) = my_Euler_thiele(t(n),V_euler_M(n),h,delta,P,S,mu);
    V_trap_M(n+1) = my_trapezoidal_thiele(t(n),V_trap_M(n),h,delta,P,S,mu);
    V_RK_M(n+1) = my_RK_thiele(t(n),V_RK_M(n),h,delta,P,S,mu);
end

% Compare reserves at year 10
% finding spot in array of year 10
YR10 = find(t==10,1);

%Reserve values by method
fprintf('Euler at year 10 (Makeham): %.10f\n', V_euler_M(YR10));
fprintf('Trapezoidal at year 10 (Makeham): %.10f\n', V_trap_M(YR10));
fprintf('RK4 at year 10 (Makeham): %.10f\n', V_RK_M(YR10));
%%
% changing mortality to Gompertz's Law
T = 20; %policy term is 20 years
delta = 0.04; %interest rate per year
P = 2500; %premiums payable continuously
S = 100000; %death benefit
mu = @(t) (2.7e-6)*(1.124).^(30+t); %mu taken in as a function - just age-dependent

% Test h = -0.05 first
N=400;
h= -T/N;
%creating equally spaced time-steps array
t = linspace(T,0,N+1); %start condition T, end condition 0, N+1 points
% t=[20, 19.95, 19.90, ... , 0]

% Setting the first term in the array as the Terminal condition = endowment benefit
V_euler_G(1) = 100000;
V_trap_G(1) = 100000;
V_RK_G(1) = 100000;

%running methods 
for n = 1:N
    V_euler_G(n+1) = my_Euler_thiele(t(n),V_euler_G(n),h,delta,P,S,mu);
    V_trap_G(n+1) = my_trapezoidal_thiele(t(n),V_trap_G(n),h,delta,P,S,mu);
    V_RK_G(n+1) = my_RK_thiele(t(n),V_RK_G(n),h,delta,P,S,mu);
end

% Compare reserves at year 10
% finding spot in array of year 10
YR10 = find(t==10,1);

%Reserve values by method
fprintf('Euler at year 10 (Gompertz): %.10f\n', V_euler_G(YR10));
fprintf('Trapezoidal at year 10 (Gompertz): %.10f\n', V_trap_G(YR10));
fprintf('RK4 at year 10 (Gompertz): %.10f\n', V_RK_G(YR10));
%%
%Compare how mortality effects within each method
fprintf('Difference in Reserves under Makeham vs Gompertz Mortality (Euler): %.10f\n', abs(V_euler_M(YR10)-V_euler_G(YR10)));
fprintf('Difference in Reserves under Makeham vs Gompertz Mortality (Trapezoidal): %.10f\n', abs(V_trap_M(YR10)-V_trap_G(YR10)));
fprintf('Difference in Reserves under Makeham vs Gompertz Mortality (RK4): %.10f\n', abs(V_RK_M(YR10)-V_RK_G(YR10)));

%% 
%test reference solution using ODE45 for non-constant mortality
%this is using MATLABs built in ODE solver 
T = 20;
S = 100000;
P = 2500;
delta = 0.04;

% Makeham mortality
mu = @(t) 0.00022 + 2.7e-6*(1.124).^(30+t);

% Thiele's ODE
f = @(t,V) delta*V + P - mu(t)*(S-V);

%testing for best tol
%options2 = odeset(RelTol=1e-10,AbsTol=1e-10);
%options3 = odeset(RelTol=1e-12,AbsTol=1e-12);
%[t_ref2,V_ref2] = ode45(f,[T 0],S,options2);
%[t_ref3,V_ref3] = ode45(f,[T 0],S,options2);
%V_reference2 = deval(sol1,10);
%V_reference3 = deval(sol2,10);
% seems to be converging as very little change between
    %  2 and 3 therefore we will use tol = 1e-12

% set up ode45. solve backwards for t=20 to t=0
options = odeset(RelTol=1e-12,AbsTol=1e-12);
[t_ref,V_ref] = ode45(f,20:-0.01:0,S,options);
% Find reference reserve at year 10
YR10_ref = find(t_ref == 10 ,1);
V_reference = V_ref(YR10_ref); 
%%%%%%%%%%%%%%%%
fprintf('Makeham reference reserve at year 10: %.10f\n',V_reference);
% Errors for h = -0.05
error_euler = abs(V_euler(YR10)-V_reference);
error_trap = abs(V_trap(YR10)-V_reference);
error_RK = abs(V_RK4(YR10)-V_reference);

% Errors for h = -0.01
error_euler01 = abs(V_euler01(YR10_2)-V_reference);
error_trap01 = abs(V_trap01(YR10_2)-V_reference);
error_RK01 = abs(V_RK01(YR10_2)-V_reference);

% Print Makeham numerical errors
fprintf('\nMakeham Mortality Error Comparison\n');
fprintf('Reference reserve at year 10: $%.10f\n',V_reference);

% Step size h = -0.05
fprintf('\nStep Size h = -0.05\n');
fprintf('Euler error: %.10f\n',error_euler);
fprintf('Trapezoidal error: %.10f\n',error_trap);
fprintf('RK4 error: %.10f\n',error_RK);

% Step size h = -0.01
fprintf('\nStep Size h = -0.01\n');
fprintf('Euler error: %.10f\n',error_euler01);
fprintf('Trapezoidal error: %.10f\n',error_trap01);
fprintf('RK4 error: %.10f\n',error_RK01);

