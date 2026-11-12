
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
% Plotting Stability bounds

t = linspace(0,T,200);

figure;
% Euler
subplot(2,1,1); %creates two stacked sub plots
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
title('Stability bound for Euler Method');
legend('Location','best');
grid on;


% RK4
subplot(2,1,2);
hold on;

for x = ages
    mu_t = B*c.^(x+t);
    hmax_RK = 2.7853./(delta + mu_t);
    plot(t,hmax_RK,'LineWidth',1.5, 'DisplayName',['Age ',num2str(x)]);

end

xlabel('Time (years)');
ylabel('Maximum |h|');
ylim([0,75]);
title('Stability Bound for RK4 Method');
legend('Location','best');
grid on;
