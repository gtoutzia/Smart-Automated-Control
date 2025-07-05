%Toutziaris Gewrgios AEM 10568
clc;
clf;
close all;

x0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
tspan = [0 20];
[t, x] = ode45(@(t, x) system_dynamics(t, x), tspan, x0);
x_actual = x(:, 1:4); % Actual system states
x_model = x(:, 9:12); % Model system states (xm)
% Plotting results
figure;

% Plot θ1 and θd1 (for first 50 seconds)
subplot(2,2,1);
plot(t, x_actual(:,1), 'r'); hold on;
plot(t, x_model(:,1), 'b--'); hold off;
title('Tracking of x1 and xm1');
xlabel('Time (s)');
ylabel('Values');
legend('x1', 'xm1');
grid on;

% Plot θ2 and θd2 (for first 50 seconds)
subplot(2,2,2);
plot(t, x_actual(:,3), 'r'); hold on;
plot(t, x_model(:,3), 'b--'); hold off;
title('Tracking of x3 and xm3');
xlabel('Time (s)');
ylabel('Values');
grid on;

% Plot error x1 - xm1
subplot(2,2,3);
plot(t, x_actual(:,1)-x_model(:,1), 'r'); hold on;
title('Tracking of error x1 - xm1');
xlabel('Time (s)');
ylabel('Values');
grid on;

% Plot error x3 - xm3
subplot(2,2,4);
plot(t, x_actual(:,3)-x_model(:,3), 'r'); hold on;
title('Tracking of error x2 - xm2');
xlabel('Time (s)');
ylabel('Values');
grid on;

figure;
% Plot θ1dot and xm2dot (for first 50 seconds)
subplot(2,2,1);
plot(t, x_actual(:,2), 'r'); hold on;
plot(t, x_model(:,2), 'b--'); hold off;
title('Tracking of x2 and xm2');
xlabel('Time (s)');
ylabel('Values');
legend('θ1dot', 'xm1dot');
grid on;

% Plot θ2dot and xm4dot (for first 50 seconds)
subplot(2,2,2);
plot(t, x_actual(:,4), 'r'); hold on;
plot(t, x_model(:,4), 'b--'); hold off;
title('Tracking of x4 and xm4');
xlabel('Time (s)');
ylabel('Values');
legend('θ2dot', 'xm2dot');
grid on;

% Plot error θ1dot - xm1dot
subplot(2,2,3);
plot(t, x_actual(:,2)-x_model(:,2), 'r'); hold on;
title('Tracking of error x2 - xm2dot');
xlabel('Time (s)');
ylabel('Values');
grid on;

% Plot error θ1dot - xm1dot
subplot(2,2,4);
plot(t, x_actual(:,4)-x_model(:,4), 'r'); hold on;
title('Tracking of error x4 - xm4dot');
xlabel('Time (s)');
ylabel('Values');
grid on;

figure;
subplot(2,1,1)
plot(t,x(:,24), 'r'); hold on;
title('u1 vs time graph')
xlabel('Time (s)');
ylabel('Values of u1');
grid on;

subplot(2,1,2)
plot(t,x(:,25), 'r'); hold on;
title('u2 vs time graph')
xlabel('Time (s)');
ylabel('Values of u2');
grid on;

function dxdt = system_dynamics(t, x)

    % Constants
    m1 = 2; m2 = 2.5;
    J1 = 0.5; J2 = 0.625;
    r = 0.5; d = 0.5; l = 0.5;
    k = 150; b = 1; g = 9.81;
    sigma0 = 1; sigma1 = 1; sigma2 = 1;
    thetasdot = 0.1; Ts = 2; Tc = 1;
    am1 = -1; am2 = -1; am3 = -1; am4 = -1;
    bm1 = 1; bm2 = 1;

    gamma = 1000; 
    sigma = 0.01;

    % Unpack state variables
    x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
    tau1 = x(5); tau2 = x(7);
    xm1 = x(9); xm2 = x(10); xm3 = x(11); xm4 = x(12);
    k1hat = x(13); k2hat = x(14); k3hat = x(15); k4hat = x(16);
    k5hat = x(17); l1hat = x(18);
    k6hat = x(19); k7hat = x(20); k8hat = x(21);
    k9hat = x(22); k10hat = x(23); l2hat = x(24);

    % Calculate theta and chi
    theta = atan2((r/2) * (cos(x3) - cos(x1)), (d + (r/2) * (sin(x1) - sin(x3))));
    chi = sqrt(d^2 + d*r*(sin(x1) - sin(x3)) + (r^2/2) * (1 - cos(x3 - x1)));
    chi_dot = (d * r * (cos(x1) * x2 - cos(x3) * x4) + (r^2 / 2) * sin(x3 - x1) * (x4 - x2)) / (2 * chi);

    % Calculate T1,T2,tau1,tau2
    T1 = sigma0 * tau1 + sigma1 * x(6) + sigma2 * x2;
    T2 = sigma0 * tau2 + sigma1 * x(8) + sigma2 * x4;

    tau1dot = x(2) - sigma0 * abs(x(2)) / (Tc + (Ts - Tc) * exp(-abs(x(2))/thetasdot)) * x(5);
    tau2dot = x(4) - sigma0 * abs(x(4)) / (Tc + (Ts - Tc) * exp(-abs(x(4))/thetasdot)) * x(7);
    
    %inputs of model
    
    r1 = (-(2*pi^3)*sin(2*pi*t)/3 - am1*pi*sin(2*pi*t)/6 - am2*pi^2*cos(2*pi*t)/3)/bm1;
    r2 = (-(pi^3)*sin(pi*t)/4 - am3*pi*sin(pi*t)/4 - am4*pi^2*cos(pi*t)/4)/bm2;  



    % Model dynamics
    xm1dot = xm2;
    xm2dot = am1 * xm1 + am2 * xm2 + r1;
    xm3dot = xm4;
    xm4dot = am3 * xm3 + am4 * xm4 + r2;

    % Errors
    e1 = x1 - xm1;
    e2 = x2 - xm2;
    e3 = x3 - xm3;
    e4 = x4 - xm4;

    %Known non-linear functions of states
    F1 = sin(x1);
    F2 = (chi - l) * r * cos(x1 - theta);
    F3 = chi_dot * r * cos(x1 - theta);

    F4 = sin(x3);
    F5 = (chi - l) * r * cos(x3 - theta);
    F6 = chi_dot * r * cos(x3 - theta);
    
    % Adaptive laws for controller parameters without σ-modification
    % (comment out to try)
%{
    k1dot = -gamma * e2 * F1;
    k2dot = -gamma * e2 * F2;
    k3dot = -gamma * e2 * F3;
    k4dot = -gamma * e2 * x1;
    k5dot = -gamma * e2 * x2;
    l1dot = -gamma * e2 * r1;

    k6dot = -gamma * e4 * F1;
    k7dot = -gamma * e4 * F2;
    k8dot = -gamma * e4 * F3;
    k9dot = -gamma * e4 * x3; 
    k10dot = -gamma * e4 * x4; 
    l2dot = -gamma * e4 * r2; 
%}
    % Adaptive laws for controller parameters (comment out to try)
   
    k1dot = -gamma * e2 * F1 -gamma*sigma*k1hat;
    k2dot = -gamma * e2 * F2 -gamma*sigma*k2hat;
    k3dot = -gamma * e2 * F3 -gamma*sigma*k3hat;
    k4dot = -gamma * e2 * x1 -gamma*sigma*k4hat;
    k5dot = -gamma * e2 * x2 -gamma*sigma*k5hat;
    l1dot = -gamma * e2 * r1 -gamma*sigma*l1hat;

    k6dot = -gamma * e4 * F4 -gamma*sigma*k6hat;
    k7dot = -gamma * e4 * F5 -gamma*sigma*k7hat;
    k8dot = -gamma * e4 * F6 -gamma*sigma*k8hat;
    k9dot = -gamma * e4 * x3 -gamma*sigma*k9hat;
    k10dot = -gamma * e4 * x4 -gamma*sigma*k10hat;
    l2dot = -gamma * e4 * r2 -gamma*sigma*l2hat;

    k1 = m1*g;
    k2 = -0.5*k;
    k3 = -0.5*b;
    k4 = m2*g;
    k5 = 0.5*k;
    k6 = 0.5*b;
    b1 = 1/J1;
    b2 = 1/J2;

    % Control inputs
    u1 = -k1hat*F1 - k2hat*F2 - k3hat*F3 + k4hat*x1 + k5hat*x2 + l1hat*r1;
    u2 = -k6hat*F4 - k7hat*F5 - k8hat*F6 + k9hat*x3 + k10hat*x4 + l2hat*r2;
    
   %{
    % System dynamics
    x1dot = x2;
    x2dot = (1/J1)*(k1hat*F1 + k2hat*F2 + k3hat*F3 - T1+ u1);
    x3dot = x4;
    x4dot = (1/J2)*(k6hat*F4 + k7hat*F5 + k8hat*F6 - T2 +u2);
   %}

     % System dynamics
    x1dot = x2;
    x2dot = b1*(k1*F1 + k2*F2 + k3*F3 - T1+ u1);
    x3dot = x4;
    x4dot = b2*(k4*F4 + k5*F5 + k6*F6 - T2 +u2);

    
    % Collect derivatives
    dxdt = [x1dot; x2dot; x3dot; x4dot; tau1dot; x(6); tau2dot; x(8); xm1dot; xm2dot; xm3dot; xm4dot; k1dot; k2dot; k3dot; k4dot; k5dot; l1dot; k6dot; k7dot; k8dot; k9dot; k10dot; l2dot;u1;u2];
end
