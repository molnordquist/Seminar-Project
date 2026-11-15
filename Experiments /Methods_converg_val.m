clear;
clc;
close all;
addpath('../SRC');
%%
%Define parameters
T = 20; 
delta = 0.04; 
P = 2500; 
S = 100000; 

mu_const = 0.01; 
mu = @(t) mu_const;

%Setting up components of Thiele's
A= delta +mu_const;
B= P-mu_const*S;

% Creating exact reserve
V_exact = @(t) -B/A + (S + B/A)*exp(A*(t-T));

%Set up changing N Values
N_vals = [10, 20, 40, 80, 160, 320, 640, 1280];

%set up error arrays to store
errors_euler = zeros(size(N_vals));
errors_RK = zeros(size(N_vals));
errors_trap = zeros(size(N_vals));

%create a for loop going through step sizes
for i = 1:length(N_vals)
    N= N_vals(i);
    h= -T/N;
    t=linspace(T,0, N+1);
    
    % initializing arrays
    V_euler = zeros(1, N+1);
    V_RK = zeros(1, N+1);
    V_trap = zeros(1, N+1);

    V_euler(1) = S;
    V_RK(1) = S;
    V_trap(1) =S;
    
    for n = 1:N
        V_euler(n+1) = my_Euler_thiele(t(n),V_euler(n),h,delta,P,S,mu);
        V_trap(n+1) = my_trapezoidal_thiele(t(n),V_trap(n),h,delta,P,S,mu);
        V_RK(n+1) = my_RK_thiele(t(n),V_RK(n),h,delta,P,S,mu);
    end
    
    yr10 = N/2 +1;
    errors_euler(i) = abs(V_euler(yr10)-V_exact(10));
    errors_RK(i) = abs(V_RK(yr10)-V_exact(10));
    errors_trap(i) = abs(V_trap(yr10)-V_exact(10));
end

%calculating order of convergence
for i= 1:length(N_vals)-1
    k_euler(i) = log(errors_euler(i)/errors_euler(i+1))/log(2);
    k_RK(i) = log(errors_RK(i)/errors_RK(i+1))/log(2);
    k_trap(i) = log(errors_trap(i)/errors_trap(i+1))/log(2);
    fprintf('Euler order of convergence: %.4f\n',k_euler);
    fprintf('RK4 order of convergence: %.4f\n',k_RK);
    fprintf('Trapezoidal order of convergence: %.4f\n',k_trap);
end

%Creatng results tables
order_converg_results = table((T./N_vals(2:end))', k_euler', k_RK', k_trap', 'VariableNames', {'h', 'Euler Order of Convergence', 'RK4 Order of Convergence', 'Trapezoidal Order of Convergence'});
disp(order_converg_results);
converg_results = table((T./N_vals)', errors_euler', errors_RK', errors_trap', 'VariableNames', {'h', 'Euler Error', 'RK4 Error', 'Trapezoidal Error'});
disp(converg_results);

%creating graphs 

%% Plot Convergence Results

h_vals = T./N_vals;

figure;
loglog(h_vals, errors_euler, 'o-', 'LineWidth', 2);
hold on;
loglog(h_vals, errors_RK, 's-', 'LineWidth', 2);
loglog(h_vals, errors_trap, '^-', 'LineWidth', 2);

xlabel('Absolute Step Size |h|');
ylabel('Absolute Error at Year 10');
title('Convergence of Numerical Methods Under Constant Mortality');
legend('Euler', 'RK4', 'Trapezoidal', 'Location', 'best');
grid on;
hold off;

%% Plot Absolute Errors

h_vals = T./N_vals;

figure;
plot(h_vals, errors_euler, 'o-', 'LineWidth', 2);
hold on;
plot(h_vals, errors_RK, 's-', 'LineWidth', 2);
plot(h_vals, errors_trap, '^-', 'LineWidth', 2);

xlabel('Absolute Step Size |h|');
ylabel('Absolute Error at Year 10');
ylim([-100 1500]);
title('Numerical Error Under Constant Mortality');
legend('Euler', 'RK4', 'Trapezoidal', 'Location', 'best');
grid on;
hold off;


