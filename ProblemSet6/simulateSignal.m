function [x_tauj, theta_tauj] = simulateSignal(tauj, Ampl, CN0, fIF, Ta,  theta0)
% simulateSignal: simulates a signal
%
% The formulation follows almost exactly page 1-2 of "notes9.pdf".
%
% x(tau_j) = A(tau_j) * D[tau_j - td(tau_j)] * cos[2*pi*fIF*tau_j + theta(tau_j)] + n(tau_j)
%
% where the only difference is that we consider the case where the
% spreading code term `C[tau_j - ts(tau_j)]` has been completely eliminated
% through "despreading" process. In other words, we assume that we found
% the code-delay that perfectly correlates a certain PRN-code sequence with
% the one in our received signal.
%
% INPUTS
%
% tauj -------- Time vector [sec]
%
% Ta ---------- Accumulation period [sec]
%
% Ampl -------- Amplitude of the signal (constant model)
%
% fIF --------- Intermediate frequency [Hz]
%
% CN0 --------- Carrier-to-Noise ratio [dB-Hz]
%
% theta0 ------ Initial phase offset [rad]
%
% OUTPUTS
%
% x    -------- simulated gnss signal as it exits the ADC under a frequency
%               plan that shifts the signal to fIF.
%
% theta ------- simulated time history of phase offsets
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+    
    % Sampling rate
    Ts = tauj(2) - tauj(1);
    
    % total number of samples
    N = length(tauj);

    % Generate amplitude of the signal at each time step
    A_tauj = simulateAmplitude(Ampl, N);

    % Generate the navigation data sequence 
    D_tauj = simulateNavigationData(N, Ts);

    % Create simulated signal.                    --- page 2 {Notes9.pdf}
    linCN0 = 10^(CN0/10);                         % linear C/N0
    varn   = (mean(A_tauj))^2/(4*linCN0*Ts);      % variance of noise
    n_tauj = sqrt(varn)*randn(1,N);               % sampled noise 

    % Generate the phase shift time history
    theta_tauj = simulatePhaseShift(tauj, theta0, 'custom');

    % Compute the signal at each tauj             --- {page 1 - Notes9.pdf}
    x_tauj = A_tauj .* D_tauj .* cos(2*pi*fIF*tauj + theta_tauj) + n_tauj; 
end

function thetaVtr = simulatePhaseShift(t, theta0, model)
% simulatePhaseShift: simulates a phase shift of the carrier 
%
%
% INPUTS
%
% t      ------------ Time vector
%
% theta0 ------------ Initial phase 
%
% model   ----------- Flag that indicates the model to implement
%
%                   Options:
%                       - 'ramp'
%                       - 'quadratic'
%
% OUTPUTS
%
% thetaVtr -------- Phase shift corresponding to each time step.
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+
    switch model

        case 'ramp'
            Ts = t(2) - t(1);
            pad_l  = round(1/Ts);
            slope  = 2*pi;
            thetaVtr  = [slope*t(1:pad_l).*ones(1,pad_l), slope*t(pad_l)*ones(1, length(t)-pad_l)];

        case 'custom'
            t1 = t(1:floor(length(t)/3));
            tF = t1(end);
            ramp1  = 2*pi*t1;
            t2 = t(floor(length(t)/3)+1: floor(2*length(t)/3)) - tF;
            tF = tF + t2(end);
            quadratic = ramp1(end) - pi*t2.^2;
            t3 = t(floor(2*length(t)/3)+1: end) - tF;
            ramp2  = quadratic(end) + 2*pi*t3; 
            
            thetaVtr = [ramp1, quadratic, ramp2];  
    end
    thetaVtr = thetaVtr + theta0;
end

function [D] = simulateNavigationData(N, Ts)
% simulateNavigationData: simulates a navigation data sequence sampled at 
%                         Ts.
%
%
% INPUTS
%
% N ------------ Number of data samples required
%
% Ts ----------- Desired sampling rate 
%
% OUTPUTS
%
%
% D -------- Navigation data 
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+
    % Rate of the navigation data 
    Td = 20e-3;
    Nd = ceil(N*Ts/Td); % number of data bits needed

    % Random sequence -1,+1. Each element of this sequence is spaced 
    % Td~20ms. Therefore, the sequence needs to be oversampled to match the 
    % sampling rate of the signal we are generating. 
    auxSeq = 2*randi([0 1], 1, Nd) - 1;

    % Oversample navigation data bits switching at 20 ms.
    D = oversampleSpreadingCode(auxSeq, Ts/Td, 0, N, Nd)';
end

function [A] = simulateAmplitude(Ampl, N)
% simulateNavigationData: simulates the amplitude of the signal
%
%
% INPUTS
%
% Ampl --------- Amplitude (constant)
%
% N ------------ Number of data bits required
%
% OUTPUTS
%
%
% A -------- Amplitude of the signal
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+
    amplitude = 40;
    A = amplitude * ones(1,N);
end