clear;
clc;
close all;
addpath('../SRC');
%%
%Gompertz Mortality

T=20; %old age so mortality behavior will increase rapidly  
delta = 0.04;
P = 2500; 
S= 100000; 

B = 2.7e-6;
c = 1.124;


%defining ages we will be looping through policy start date ages
ages = [30 50 70 90];
%as defiing N values to change step-sizes
N_vals = [1,2,4,10,20];

%define matrix to save values 
Gompertz_euler = zeros(length(ages),length(N_vals));
Gompertz_RK = zeros(length(ages),length(N_vals));
Gompertz_trap = zeros(length(ages),length(N_vals));
Gompertz_exact = zeros(size(ages));

%loop through policy starting age
colors = lines(length(ages));
figure;
hold on;
for j = 1:length(ages)
    x = ages(j);
    %mortality must be in this loop because it changes with age
    mu = @(t) B*c.^(x+t);

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
    Gompertz_euler(j,:) = V0_euler;
    Gompertz_RK(j,:) = V0_RK;
    Gompertz_trap(j,:) = V0_trap;

    %calculate refrence solution for Gompertz
    I = @(t) exp(-delta*t - (B*c.^(x+t))/log(c));
    full_int_gomp = @(s) I(s).*(P-S*mu(s));
    integral_sol_gomp = integral(full_int_gomp,0,T);
    Gompertz_exact(j) = (I(T)*S - integral_sol_gomp)/I(0);

    % Plot numerical reserves for this age
    plot(T./N_vals,V0_euler,'o-', ...
        'Color',colors(j,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('Euler at Age %d',x));

    plot(T./N_vals,V0_RK,'s-', ...
        'Color',colors(j,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('RK4 at Age %d',x));

    plot(T./N_vals,V0_trap,'^-', ...
        'Color',colors(j,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('Trapezoidal at Age %d',x));
    
    yline(Gompertz_exact(j), '--', ...
    'Color',colors(j,:), ...
    'LineWidth',2, ...
    'DisplayName',sprintf('Exact Solution Expected at age %d',ages(j)));
end

xlabel('Absolute Step Size |h|');
ylabel('Reserve at t = 0 ($)');
ylim([-30000 350000])
title('Numerical Reserves by Age under Gompertz Mortality');
legend('Location','eastoutside');
grid on;
    
    %results = table(N_vals', (T./N_vals)',V0_euler', V0_RK', V0_trap', ...
        %'VariableNames', {'N','StepSize', ...
        %'EulerReserve','RK4Reserve','TrapezoidalReserve'});
    
    %disp(results);
    
   
    

% Find reserve at year 10, indicated year of textbook. Same as 40 years old
%YR0 = find(t == 0,1);
%fprintf('Gompertz reserve at year 0: %.10f\n',Gompertz_exact(YR0));

%error_euler = abs(V0_euler - Gompertz_exact(YR0));
%error_RK = abs(V0_RK - Gompertz_exact(YR0));
%error_trap = abs(V0_trap - Gompertz_exact(YR0));

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



%% Makeham 
T = 20;
S = 100000;
P = 2500;
delta = 0.04;
x = 70;

% Seting up Makeham mortality parameters
A = 0.00022;
B = 2.7e-6;
c = 1.124;

%defining ages we will be looping through policy start date ages
ages = [30 50 70 90];
%as defiing N values to change step-sizes
N_vals = [1,2,4,10,20];
%initiatling 
Makeham_euler = zeros(length(ages),length(N_vals));
Makeham_RK = zeros(length(ages),length(N_vals));
Makeham_trap = zeros(length(ages),length(N_vals));
Makeham_exact = zeros(size(ages));

colors = lines(length(ages));
figure;
hold on;
%loop through policy starting age
for j = 1:length(ages)
    x = ages(j);
    %mortality must be in this loop because it changes with age
    mu = @(t) A + B*c.^(x+t);

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
        V0_euler_m(i) = V_euler(end);
        V0_RK(i) = V_RK(end);
        V0_trap(i) = V_trap(end);

    end
    %updates the initiated arrays
    Makeham_euler(j,:) = V0_euler;
    Makeham_RK(j,:) = V0_RK;
    Makeham_trap(j,:) = V0_trap;

    % Calculating the refrence solution for Makeham mortality
    I = @(t) exp(-delta*t -A*t - (B*c.^(x+t))/log(c));
    full_int_mak = @(s) I(s).*(P-S*mu(s));
    integral_sol_mak = integral(full_int_mak,0,T);
    Makeham_exact(j) = (I(T)*S - integral_sol_mak)/I(0);

    plot(T./N_vals,V0_euler,'o-', ...
        'Color',colors(j,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('Euler at Age %d',x));

    plot(T./N_vals,V0_RK,'s-', ...
        'Color',colors(j,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('RK4 at Age %d',x));

    plot(T./N_vals,V0_trap,'^-', ...
        'Color',colors(j,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('Trapezoidal at Age %d',x));
    yline(Makeham_exact(j), '--', ...
    'Color',colors(j,:), ...
    'LineWidth',2, ...
    'DisplayName',sprintf('Exact Solution Expected at age %d',ages(j)));
end

xlabel('Absolute Step Size |h|');
ylabel('Reserve at t = 0 ($)');
ylim([-30000 350000])
title('Numerical Reserves by Age under Makeham Mortality');
legend('Location','eastoutside');
grid on;




% Find reserve at year 10, indicated year of textbook. Same as 40 years old
%YR0 = find(t == 0,1);
%fprintf('Makeham reserve at year 0: %.10f\n',Makeham_exact(YR0));

%error_euler = abs(V0_euler - Makeham_exact(YR0));
%error_RK = abs(V0_RK - Makeham_exact(YR0));
%error_trap = abs(V0_trap - Makeham_exact(YR0));

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
%YR0 = find(t == 0,1);
%fprintf('Constant reserve at year 0: %.10f\n',V_constant);

%error_euler = abs(V0_euler - V_constant);
%error_RK = abs(V0_RK - V_constant);
%error_trap = abs(V0_trap - V_constant);

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

% Save results to plot together
Constant_euler = V0_euler;
Constant_RK = V0_RK;
Constant_trap = V0_trap;
Constant_exact = V_constant;

%%
%Plotting Absolute error

figure;

% Gompertz mortality

h_vals = T./N_vals;

for j = 1:length(ages)

    subplot(2,4,j);
    hold on;

    errors_euler = abs(Gompertz_euler(j,:)-Gompertz_exact(j));
    errors_RK = abs(Gompertz_RK(j,:)-Gompertz_exact(j));
    errors_trap = abs(Gompertz_trap(j,:)-Gompertz_exact(j));
    
    plot(h_vals,errors_euler,'Color', 'r','LineWidth', 1.5);
    plot(h_vals,errors_RK,'Color', 'b','LineWidth', 1.5);
    plot(h_vals,errors_trap,'LineWidth', 1.5);

    xlabel('Step Size |h|');
    ylabel('Absolute Error (Gompertz Mortality)');
    title(['Age ',num2str(ages(j)), ' (Gompertz)']);

    legend('Euler','RK4','Trapezoidal');
    grid on;
%calculating order of convergence
    for i= 1:length(N_vals)-1
        k_euler(i) = log(errors_euler(i)/errors_euler(i+1))/log(h_vals(i)/h_vals(i+1));
        k_RK(i) = log(errors_RK(i)/errors_RK(i+1))/log(h_vals(i)/h_vals(i+1));
        k_trap(i) = log(errors_trap(i)/errors_trap(i+1))/log(h_vals(i)/h_vals(i+1));
    end

%Creatng results tables
    order_converg_results = table((h_vals(2:end))', k_euler', k_RK', k_trap', 'VariableNames', {'h', 'Euler Order of Convergence', 'RK4 Order of Convergence', 'Trapezoidal Order of Convergence'});
    disp(order_converg_results);
    converg_results = table(h_vals', errors_euler', errors_RK', errors_trap', 'VariableNames', {'h', 'Euler Error', 'RK4 Error', 'Trapezoidal Error'});
    disp(converg_results);
end



% Makeham mortality
for j = 1:length(ages)
    
    subplot(2,4,j+4);
    hold on;
    
    errors_euler = abs(Makeham_euler(j,:)-Makeham_exact(j));
    errors_RK = abs(Makeham_RK(j,:)-Makeham_exact(j));
    errors_trap = abs(Makeham_trap(j,:)-Makeham_exact(j));

    plot(h_vals,errors_euler,'Color', 'r','LineWidth', 1.5);
    plot(h_vals,errors_RK,'Color', 'b', 'LineWidth', 1.5);
    plot(h_vals,errors_trap, 'LineWidth', 1.5);
    xlabel('Step Size |h|');
    ylabel('Absolute Error');
    title(['Age ',num2str(ages(j)), ' (Makeham)']);

    legend('Euler','RK4','Trapezoidal');
    grid on;
%calculating order of convergence
    for i= 1:length(N_vals)-1
        k_euler(i) = log(errors_euler(i)/errors_euler(i+1))/log(h_vals(i)/h_vals(i+1));
        k_RK(i) = log(errors_RK(i)/errors_RK(i+1))/log(h_vals(i)/h_vals(i+1));
        k_trap(i) = log(errors_trap(i)/errors_trap(i+1))/log(h_vals(i)/h_vals(i+1));
    end

%Creatng results tables
    order_converg_results = table((h_vals(2:end))', k_euler', k_RK', k_trap', 'VariableNames', {'h', 'Euler Order of Convergence', 'RK4 Order of Convergence', 'Trapezoidal Order of Convergence'});
    disp(order_converg_results);
    converg_results = table(h_vals', errors_euler', errors_RK', errors_trap', 'VariableNames', {'h', 'Euler Error', 'RK4 Error', 'Trapezoidal Error'});
    disp(converg_results);
end

sgtitle('Absolute Error vs. Step Size');




