%%Creating function files 
% Experiment 1: Comparing Euler, Trapezoidal, and RK4
clear;
clc;
close all;

% Make sure MATLAB can find the functions in SRC
addpath('../SRC');

%% DefineParameters
T = 10; %terminal condition
delta = 0.03;
P = 800;
S = 100000;

%% Define constant mortality function
mu = @(t) 0.01;
% Number of steps
N = 100;
% Negative step size (working backward becuase we define eulers to take backwards steps)
h = -T/N;
% Time grid
t = linspace(T,0,N+1);

% Initialize numerical solutions
V_euler = zeros(1,N+1);
V_trap = zeros(1,N+1);
V_RK = zeros(1,N+1);

% Terminal conditions: V(T) = 0
V_euler(1) = 0;
V_trap(1) = 0;
V_RK(1) = 0;

%% Run numerical methods
for n = 1:N
    V_euler(n+1) = my_Euler_thiele(t(n),V_euler(n),h,delta,P,S,mu);

    V_trap(n+1) = my_trapezoidal_thiele(t(n),V_trap(n),h,delta,P,S,mu);

    V_RK(n+1) = my_RK_thiele(t(n),V_RK(n),h,delta,P,S,mu);
end

%% Plot
figure;
plot(t,V_euler,'LineWidth',1.5);
hold on;

plot(t,V_trap,'--','LineWidth',1.5);
plot(t,V_RK,':','LineWidth',2);

xlabel('Time (years)');
ylabel('Reserve ($)');
title('Comparison of Numerical Methods for Thiele''s ODE');

legend('Euler','Trapezoidal','RK4');
grid on;