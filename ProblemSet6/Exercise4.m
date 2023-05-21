clc; close all; clear all;

Bn = 10;

t=(0:0.01:10)';
step_u = ones(size(t));
ramp_u = t;
quad_u = 0.5*t.^2;

%% First order
% Tune such that the loop noise bandwidth is Bn
K = 4*Bn;
% Generate the system
num = K;
den = [1, K];
first_order = tf(num,den);

[y_step, t_step] = lsim(first_order, step_u, t);
[y_ramp, t_ramp] = lsim(first_order, ramp_u, t);
[y_quad, t_quad] = lsim(first_order, quad_u, t);

figure()
subplot(1,3,1)
    plot(t, step_u - y_step, LineWidth=2)
    title('Step input error of the first order system')
    xlabel('time [s]')
subplot(1,3,2)
    plot(t, ramp_u - y_ramp, LineWidth=2)
    title('Ramp input error of the first order system')
    xlabel('time [s]')
subplot(1,3,3)
    plot(t, quad_u - y_quad, LineWidth=2)
    title('Quadratic input error of the first order system')
    xlabel('time [s]')

% Frequency response
figure()
bode(first_order)
title('Frequency response of the first order system')

%% Second order
K = (8/3)*Bn;
a = K/2;
% Generate the system
num = K*[1, a];
den = [1, K, K*a];
second_order = tf(num,den);

[y_step, t_step] = lsim(second_order, step_u, t);
[y_ramp, t_ramp] = lsim(second_order, ramp_u, t);
[y_quad, t_quad] = lsim(second_order, quad_u, t);

figure()
subplot(1,3,1)
    plot(t, step_u - y_step, LineWidth=2)
    title('Step input error of the second order system')
    xlabel('time [s]')
subplot(1,3,2)
    plot(t, ramp_u - y_ramp, LineWidth=2)
    title('Ramp input error of the second order system')
    xlabel('time [s]')
subplot(1,3,3)
    plot(t, quad_u - y_quad, LineWidth=2)
    title('Quadratic input error of the second order system')
    xlabel('time [s]')

% Frequency response
figure()
bode(second_order)
title('Frequency response of the second order system')

%% Third order
a = 1.2*Bn;
b = a^2/2;
K = 2*a;
% Generate the system
num = K*[1, a, b];
den = [1, K, K*a, K*b];
third_order = tf(num,den);

[y_step, t_step] = lsim(third_order, step_u, t);
[y_ramp, t_ramp] = lsim(third_order, ramp_u, t);
[y_quad, t_quad] = lsim(third_order, quad_u, t);

figure()
subplot(1,3,1)
    plot(t, step_u - y_step, LineWidth=2)
    title('Step input error of the third order system')
    xlabel('time [s]')
subplot(1,3,2)
    plot(t, ramp_u - y_ramp, LineWidth=2)
    title('Ramp input error of the third order system')
    xlabel('time [s]')
subplot(1,3,3)
    plot(t, quad_u - y_quad, LineWidth=2)
    title('Quadratic input error of the third order system')
    xlabel('time [s]')

% Frequency response
figure()
bode(third_order)
title('Frequency response of the third order system')