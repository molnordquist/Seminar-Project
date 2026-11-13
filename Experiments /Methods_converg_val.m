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
YR0 = find(t==0,1);
YR10 = find(t==10,1);
%using coded constant exact solution
V_exact_0 = V_exact_constant(0);
V_exact_10 = V_exact_constant(10);
error_euler_0 = abs(V_euler(YR10)-V_exact_0);
error_RK_0 = abs(V_RK(YR0)-V_exact_0);
error_trap_0 = abs(V_trap(YR0)-V_exact_0);
error_euler = abs(V_euler(YR0)-V_exact_10);
error_RK = abs(V_RK(YR10)-V_exact_10);
error_trap = abs(V_trap(YR10)-V_exact_10);

%Reserve values by method
fprintf('Step size = %.4f\n', h)
fprintf('Error of Euler at year 10: %.10f\n', error_euler);
fprintf('Error of RK4 at year 10: %.10f\n', error_RK);
fprintf('Error of Trapezoidal at year 10: %.10f\n', error_trap);
fprintf('Error of Euler at year 0: %.10f\n', error_euler_0);
fprintf('Error of RK4 at year 0: %.10f\n', error_RK_0);
fprintf('Error of Trapezoidal at year 0: %.10f\n', error_trap_0);

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
% New step-size -0.0125

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
N_vals = [10 20, 40, 80, 160, 320, 640, 1280];

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

