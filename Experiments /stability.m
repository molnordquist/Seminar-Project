clear;
clc;
close all;
addpath('../SRC');
%% Stability Tests 

%Parameters
delta = 0.04;
T = 20;
x = 30;

% Mortality parameters
mu_const = 0.01;
A = 0.00022;
B = 2.7e-6;
c = 1.124;

%%
% Constant mortality 

%stability conditions
hmax_euler = 2/abs(delta + mu_const);
hmax_RK = 2.7853/abs(delta + mu_const);

fprintf('Euler (Constant) at age 30: |h| < %.10f\n',hmax_euler);
fprintf('RK4 (Constant) at age 30: |h| < %.10f\n',hmax_RK);

%% 
% Makeham Mortality model
mu_T = A + B*c^(x+T);

%stability conditions
hmax_euler = 2/abs(delta + mu_T);
hmax_RK = 2.7853/abs(delta + mu_T);

fprintf('Euler (Makeham) at age 30: |h| < %.10f\n',hmax_euler);
fprintf('RK4 (Makeham)at age 30: |h| < %.10f\n',hmax_RK);

%% 
% Gompertz Mortality model
mu_T = B*c^(x+T);

%stability conditions
hmax_euler = 2/abs(delta + mu_T);
hmax_RK = 2.7853/abs(delta + mu_T);

fprintf('Euler (Gompertz) at age 30: |h| < %.10f\n',hmax_euler);
fprintf('RK4 (Gompertz) at age 30: |h| < %.10f\n',hmax_RK);

%% 
% Stability for Different Ages of policy start date (x)
% Gompertz Mortality

ages = [20 30 40 50 60 70 80 90];

for x = ages
    mu_T = B*c^(x+T);

    Gomp_hmax_euler = 2/abs(delta + mu_T);
    Gomp_hmax_RK = 2.7853/abs(delta + mu_T);

    fprintf('Age %d:\n',x);
    fprintf('Euler |h| < %.4f\n',Gomp_hmax_euler);
    fprintf('RK4   |h| < %.4f\n\n',Gomp_hmax_RK);

end

%% 
% Plotting Stability bounds under all mortality

t = linspace(0,T,200);

figure('Position',[100,100,1100,1000]);

% Euler for constant mortality 
subplot(3,2,1); %creates stacked subplots
mu_t = mu_const;
%calculating all stability bounds for all mortality vals
hmax_euler = 2./(delta + mu_t); 
yline(hmax_euler, 'LineWidth', 1.5)
    
xlabel('Time (years)');
ylabel('Maximum |h|');
ylim([0,75]);
title('Stability bound for Euler Method under Constant mortality');
legend('Location','best');
grid on;

% RK4 for constant mortality
subplot(3,2,2);
mu_t = mu_const;
hmax_RK = 2.7853./(delta + mu_t);
yline(hmax_RK, 'LineWidth', 1.5)

xlabel('Time (years)');
ylabel('Maximum |h|');
ylim([0,75]);
title('Stability Bound for RK4 Method under Constant Mortality');
legend('Location','best');
grid on;

%%%%
%Euler for Gompertz Mortality
subplot(3,2,3)
hold on;
for x = ages
    %calculates mortality for age x for all time steps
    mu_t = B*c.^(x+t);
    %calculating all stability bounds for all mortality vals
    hmax_euler = 2./(delta + mu_t); 
    %num2str concerts numbers to a character array
    plot(t,hmax_euler,'LineWidth',1.5,'DisplayName',['Age ',num2str(x)]);
end

xlabel('Time (years)');
ylabel('Maximum |h|');
ylim([0,75]);
title('Stability bound for Euler Method under Gompertz mortality');
legend('Location','best');
grid on;

% RK4 for Gompertz Mortality
subplot(3,2,4);
hold on;

for x = ages
    mu_t = B*c.^(x+t);
    hmax_RK = 2.7853./(delta + mu_t);
    plot(t,hmax_RK,'LineWidth',1.5, 'DisplayName',['Age ',num2str(x)]);

end

xlabel('Time (years)');
ylabel('Maximum |h|');
ylim([0,75]);
title('Stability Bound for RK4 Method under Gompertz Mortality');
legend('Location','best');
grid on;

%%%%
% Euler for Makeham Mortality
subplot(3,2,5); %creates stacked sub plots
hold on;
for x = ages
    %calculates mortality for age x for all time steps
    mu_t = A+ B*c.^(x+t);
    %calculating all stability bounds for all mortality vals
    hmax_euler = 2./(delta + mu_t); 
    %num2str concerts numbers to a character array
    plot(t,hmax_euler,'LineWidth',1.5,'DisplayName',['Age ',num2str(x)]);
end

xlabel('Time (years)');
ylabel('Maximum |h|');
ylim([0,75]);
title('Stability bound for Euler Method under Makeham Mortality');
legend('Location','best');
grid on;


% RK4 for Makeham Mortality
subplot(3,2,6);
hold on;

for x = ages
    mu_t = A+ B*c.^(x+t);
    hmax_RK = 2.7853./(delta + mu_t);
    plot(t,hmax_RK,'LineWidth',1.5, 'DisplayName',['Age ',num2str(x)]);

end

xlabel('Time (years)');
ylabel('Maximum |h|');
ylim([0,75]);
title('Stability Bound for RK4 Method under Makeham Mortality');
legend('Location','best');
grid on;

%%
%Gompertz Mortality

T=20;
x= 70; %old age so mortality behavior will increase rapidly  
delta = 0.04;
P = 2500; 
S= 100000; 

B = 2.7e-6;
c = 1.124;
mu = @(t) B*c^(x+T);
mu_T = mu(T);

N_vals = [1,2,4,10,20];

V0_euler = zeros(size(N_vals));
V0_RK = zeros(size(N_vals));
V0_trap = zeros(size(N_vals));

for i = 1:length(N_vals)
    N = N_vals(i);
    h = -T/N;
    t = linspace(T,0,N+1);

    % Initialize the reserves
    V_euler = zeros(1,N+1);
    V_RK = zeros(1,N+1);
    V_trap = zeros(1,N+1);

    %creating first entry of array as terminal condition
    V_euler(1) = S;
    V_RK(1) = S;
    V_trap(1) = S;

    % Solve using numerical methods
    for n = 1:N
        V_euler(n+1) = my_Euler_thiele(t(n),V_euler(n),h,delta,P,S,mu);
        V_RK(n+1) = my_RK_thiele(t(n),V_RK(n),h,delta,P,S,mu);
        V_trap(n+1) = my_trapezoidal_thiele(t(n),V_trap(n),h,delta,P,S,mu);
    end

    % Store reserve at t = 0 for each step size
    V0_euler(i) = V_euler(end);
    V0_RK(i) = V_RK(end);
    V0_trap(i) = V_trap(end);

end


results = table(N_vals', (T./N_vals)',V0_euler', V0_RK', V0_trap', ...
    'VariableNames', {'N','StepSize', ...
    'EulerReserve','RK4Reserve','TrapezoidalReserve'});

disp(results);

%recalculating exact under gompertz. Creating derived Integrating factor
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
fprintf('Gompertz reserve at year 0: %.10f\n',V_exact_gompertz(YR0));

error_euler = abs(V0_euler - V_exact_gompertz(YR0));
error_RK = abs(V0_RK - V_exact_gompertz(YR0));
error_trap = abs(V0_trap - V_exact_gompertz(YR0));

%figure;

%plot(T./N_vals,V0_euler,'o-','LineWidth',1.5);
%hold on;
%plot(T./N_vals,V0_RK,'s-','LineWidth',1.5);
%plot(T./N_vals,V0_trap,'^-','LineWidth',1.5);

%xlabel('Absolute Step Size |h|');
%ylabel('Reserve at t = 0 ($)');
%title('Numerical Reserves under Gompertz Mortality');
%legend('Euler','RK4','Trapezoidal','Location','best');
%grid on;

%Save results to plot together 
Gompertz_euler = V0_euler;
Gompertz_RK = V0_RK; 
Gompertz_trap = V0_trap;
Gompertz_exact = V_exact_gompertz(1);

%% Makeham 
T = 20;
S = 100000;
P = 2500;
delta = 0.04;
x = 70;

% Seting up Gompertz mortality parameters
A = 0.00022;
B = 2.7e-6;
c = 1.124;
mu = @(t) A + B*c.^(x+t);
N_vals = [1,2,4,10,20];

V0_euler = zeros(size(N_vals));
V0_RK = zeros(size(N_vals));
V0_trap = zeros(size(N_vals));

for i = 1:length(N_vals)
    N = N_vals(i);
    h = -T/N;
    t = linspace(T,0,N+1);

    % Initialize the reserves
    V_euler = zeros(1,N+1);
    V_RK = zeros(1,N+1);
    V_trap = zeros(1,N+1);

    %creating first entry of array as terminal condition
    V_euler(1) = S;
    V_RK(1) = S;
    V_trap(1) = S;

    % Solve using numerical methods
    for n = 1:N
        V_euler(n+1) = my_Euler_thiele(t(n),V_euler(n),h,delta,P,S,mu);
        V_RK(n+1) = my_RK_thiele(t(n),V_RK(n),h,delta,P,S,mu);
        V_trap(n+1) = my_trapezoidal_thiele(t(n),V_trap(n),h,delta,P,S,mu);
    end

    % Store reserve at t = 0 for each step size
    V0_euler(i) = V_euler(end);
    V0_RK(i) = V_RK(end);
    V0_trap(i) = V_trap(end);

end


results = table(N_vals', (T./N_vals)',V0_euler', V0_RK', V0_trap', ...
    'VariableNames', {'N','StepSize', ...
    'EulerReserve','RK4Reserve','TrapezoidalReserve'});

disp(results);
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
fprintf('Makeham reserve at year 0: %.10f\n',V_exact_makeham(YR0));

error_euler = abs(V0_euler - V_exact_makeham(YR0));
error_RK = abs(V0_RK - V_exact_makeham(YR0));
error_trap = abs(V0_trap - V_exact_makeham(YR0));

%figure;

%plot(T./N_vals,V0_euler,'o-','LineWidth',1.5);
%hold on;
%plot(T./N_vals,V0_RK,'s-','LineWidth',1.5);
%plot(T./N_vals,V0_trap,'^-','LineWidth',1.5);

%xlabel('Absolute Step Size |h|');
%ylabel('Reserve at t = 0 ($)');
%title('Numerical Reserves under Makeham Mortality');
%legend('Euler','RK4','Trapezoidal','Location','best');
%grid on;

%Save results to plot together 
Makeham_euler = V0_euler;
Makeham_RK = V0_RK; 
Makeham_trap = V0_trap;
Makeham_exact = V_exact_makeham(1);


%%
% Constant Mortality 
T = 20;
S = 100000;
P = 2500;
delta = 0.04;
x = 70; %age

% mortality
mu_const = 0.01;
mu = @(t) mu_const;

%creating timegrid
N_vals = [1,2,4,10,20];

V0_euler = zeros(size(N_vals));
V0_RK = zeros(size(N_vals));
V0_trap = zeros(size(N_vals));

for i = 1:length(N_vals)
    N = N_vals(i);
    h = -T/N;
    t = linspace(T,0,N+1);
    % Initialize the reserves
    V_euler = zeros(1,N+1);
    V_RK = zeros(1,N+1);
    V_trap = zeros(1,N+1);

    %creating first entry of array as terminal condition
    V_euler(1) = S;
    V_RK(1) = S;
    V_trap(1) = S;

    % Solve using numerical methods
    for n = 1:N
        V_euler(n+1) = my_Euler_thiele(t(n),V_euler(n),h,delta,P,S,mu);
        V_RK(n+1) = my_RK_thiele(t(n),V_RK(n),h,delta,P,S,mu);
        V_trap(n+1) = my_trapezoidal_thiele(t(n),V_trap(n),h,delta,P,S,mu);
    end

    % Store reserve at t = 0 for each step size
    V0_euler(i) = V_euler(end);
    V0_RK(i) = V_RK(end);
    V0_trap(i) = V_trap(end);

end


results = table(N_vals', (T./N_vals)',V0_euler', V0_RK', V0_trap', ...
    'VariableNames', {'N','StepSize', ...
    'EulerReserve','RK4Reserve','TrapezoidalReserve'});

disp(results);
%Setting up components
A= delta +mu_const;
B= P-mu_const*S;

% Creating exact reserve
V_exact_constant = @(t) -B/A + (S + B/A)*exp(A*(t-T));
V_constant = V_exact_constant(0);

% Find reserve at year 10, indicated year of textbook. Same as 40 years old
YR0 = find(t == 0,1);
fprintf('Constant reserve at year 0: %.10f\n',V_constant);

error_euler = abs(V0_euler - V_constant);
error_RK = abs(V0_RK - V_constant);
error_trap = abs(V0_trap - V_constant);

%figure;

%plot(T./N_vals,V0_euler,'o-','LineWidth',1.5);
%hold on;
%plot(T./N_vals,V0_RK,'s-','LineWidth',1.5);
%plot(T./N_vals,V0_trap,'^-','LineWidth',1.5);
%yline(V_constant, 'k--', 'LineWidth', 1.5);
%xlabel('Absolute Step Size |h|');
%ylabel('Reserve at t = 0 ($)');
%title('Numerical Reserves under Constant Mortality');
%legend('Euler','RK4','Trapezoidal','Location','best');

%Save results to plot together 
Constant_euler = V0_euler;
Constant_RK = V0_RK; 
Constant_trap = V0_trap;
Constant_exact = V_constant(1);

%%
%Plotting all three mortality models together

h_vals = T./N_vals;

figure('Position',[100,100,1500,450]);

% Constant Mortality - Left
subplot(1,3,1);
plot(h_vals,Constant_euler,'o-','LineWidth',1.5);
hold on;
plot(h_vals,Constant_RK,'s-','LineWidth',1.5);
plot(h_vals,Constant_trap,'^-','LineWidth',1.5);
yline(Constant_exact,'k--','LineWidth',1.5);

xlabel('Absolute Step Size |h|');
ylabel('Reserve at t = 0 ($)');
title('Constant Mortality');
legend('Euler','RK4','Trapezoidal','Reference','Location','best');
grid on;

% Gompertz Mortality - Middle
subplot(1,3,2);
plot(h_vals,Gompertz_euler,'o-','LineWidth',1.5);
hold on;
plot(h_vals,Gompertz_RK,'s-','LineWidth',1.5);
plot(h_vals,Gompertz_trap,'^-','LineWidth',1.5);
yline(Gompertz_exact,'k--','LineWidth',1.5);

xlabel('Absolute Step Size |h|');
ylabel('Reserve at t = 0 ($)');
title('Gompertz Mortality');
legend('Euler','RK4','Trapezoidal','Reference','Location','best');
grid on;

% Makeham Mortality - Right
subplot(1,3,3);
plot(h_vals,Makeham_euler,'o-','LineWidth',1.5);
hold on;
plot(h_vals,Makeham_RK,'s-','LineWidth',1.5);
plot(h_vals,Makeham_trap,'^-','LineWidth',1.5);
yline(Makeham_exact,'k--','LineWidth',1.5);

xlabel('Absolute Step Size |h|');
ylabel('Reserve at t = 0 ($)');
title('Makeham Mortality');
legend('Euler','RK4','Trapezoidal','Reference','Location','best');
grid on;



