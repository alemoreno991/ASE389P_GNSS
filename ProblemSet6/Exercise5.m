clc; close all; clear all;

Bn = 10;

% Number controlled oscillator
NCOs = tf(1,[1,0]);

% First order
% Tune such that the loop noise bandwidth is Bn
K = 4*Bn;
num = K;
den = 1;
Ds_1 = tf(num,den);


% Second order
K = (8/3)*Bn;
a = K/2;
% Generate the system
num = K*[1, a];
den = [1, 0];
Ds_2 = tf(num,den);

% Third order
a = 1.2*Bn;
b = a^2/2;
K = 2*a;
% Generate the system
num = K*[1, a, b];
den = [1, 0, 0];
Ds_3 = tf(num,den);


%% Discretize
T = [1, 10, 20, 40]*1e-3; % Discretization period

for ii = 1:length(T)
    BnT = Bn*T(ii)  % Eyeballing: 0.1 works relatively fine 
                    % (smaller even better!)

    % Phase averaging transfer function
    Az = tf([1 1],[2,0], T(ii));
    
    % Discretized number-controlled oscillator
    NCOz = c2d(NCOs,T(ii),'zoh');
    
    % First order
    Dz_1 = c2d(Ds_1, T(ii), 'zoh');
    Hz_1{ii} = feedback(Az*Dz_1*NCOz, 1);
    
    
    % Second order
    Dz_2 = c2d(Ds_2, T(ii), 'zoh');
    Hz_2{ii} = feedback(Az*Dz_2*NCOz, 1);

    % Third order
    Dz_3 = c2d(Ds_3, T(ii), 'zoh');
    Hz_3{ii} = feedback(Az*Dz_3*NCOz, 1);
end

%% PLOTS

for ii = 1:length(T)
    t{ii} = (0:T(ii):1)';
    step_u{ii} = ones(size(t{ii}));
    ramp_u{ii} = t{ii};
    quad_u{ii} = t{ii}.^2;
end

% First order system
figure()    
    subplot(1,3,1)
        legend_str = [];
        for ii = 1:length(T)
            legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
            lsimplot(Hz_1{ii}, step_u{ii}, t{ii})
            hold on
        end
        hold off
        legend(legend_str)
        title('Step input error of the first order system')
        xlabel('time [s]')
    subplot(1,3,2)
        legend_str = [];
        for ii = 1:length(T)
            legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
            lsimplot(Hz_1{ii}, ramp_u{ii}, t{ii})
            hold on
        end
        hold off
        legend(legend_str)
        title('Ramp input error of the first order system')
        xlabel('time [s]')
    subplot(1,3,3)
        legend_str = [];
        for ii = 1:length(T)
            legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
            lsimplot(Hz_1{ii}, quad_u{ii}, t{ii})
            hold on
        end
        hold off
        legend(legend_str)
        title('Quadratic input error of the first order system')
        xlabel('time [s]')

% Second order system
figure()    
    subplot(1,3,1)
        legend_str = [];
        for ii = 1:length(T)
            legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
            lsimplot(Hz_2{ii}, step_u{ii}, t{ii})
            hold on
        end
        hold off
        legend(legend_str)
        title('Step input error of the second order system')
        xlabel('time [s]')
    subplot(1,3,2)
        legend_str = [];
        for ii = 1:length(T)
            legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
            lsimplot(Hz_2{ii}, ramp_u{ii}, t{ii})
            hold on
        end
        hold off
        legend(legend_str)
        title('Ramp input error of the second order system')
        xlabel('time [s]')
    subplot(1,3,3)
        legend_str = [];
        for ii = 1:length(T)
            legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
            lsimplot(Hz_2{ii}, quad_u{ii}, t{ii})
            hold on
        end
        hold off
        legend(legend_str)
        title('Quadratic input error of the second order system')
        xlabel('time [s]')

% Third order system
figure()    
    subplot(1,3,1)
        legend_str = [];
        for ii = 1:length(T)
            legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
            lsimplot(Hz_3{ii}, step_u{ii}, t{ii})
            hold on
        end
        hold off
        legend(legend_str)
        title('Step input error of the third order system')
        xlabel('time [s]')
    subplot(1,3,2)
        legend_str = [];
        for ii = 1:length(T)
            legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
            lsimplot(Hz_3{ii}, ramp_u{ii}, t{ii})
            hold on
        end
        hold off
        legend(legend_str)
        title('Ramp input error of the third order system')
        xlabel('time [s]')
    subplot(1,3,3)
        legend_str = [];
        for ii = 1:length(T)
            legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
            lsimplot(Hz_3{ii}, quad_u{ii}, t{ii})
            hold on
        end
        hold off
        legend(legend_str)
        title('Quadratic input error of the third order system')
        xlabel('time [s]')        
%% Frequency response
ii = 4;

% Calculate Bn_act
%--------------------------------------------------------------------------
% This basically numerically computes the expression of Bn
% 
% Bn = (1/2*pi) * int_{-\infty}^{\infty} ( |H(j*w)|^2 dw )
%
% where we assumed (H(0) = 1).
%
% Therefore, numerically
%
% Bn_act = (1/2*pi) * sum( |H(j*w)|^2 dw )
%
% Finally, due to the assumption of H(0) = 1 we normalize to make sure the 
% expression holds even for H(0) ~= 1.
%--------------------------------------------------------------------------
walias = pi/T(ii);
wvec = [0:10000]'*(walias/10000);
[magvec,~] = bode(Hz_1{ii}, wvec);
magvec = magvec(:);
Bn_act = sum(magvec.^2)*mean(diff(wvec))/(2*pi*(magvec(1,1)^2))


%% Bode plots

% First order 
figure()
legend_str = [];
for ii = 1:length(T)
    legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
    bode(Hz_1{ii})
    hold on
end
hold off
title('Frequency response of the first order system')
legend(legend_str)

% Second order 
figure()
legend_str = [];
for ii = 1:length(T)
    legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
    bode(Hz_2{ii})
    hold on
end
hold off
title('Frequency response of the second order system')
legend(legend_str)

% Third order 
figure()
legend_str = [];
for ii = 1:length(T)
    legend_str = [legend_str, "T = " + string(T(ii)*1000) + " [ms]"];
    bode(Hz_3{ii})
    hold on
end
hold off
title('Frequency response of the third order system')
legend(legend_str)