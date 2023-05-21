clc; close all;clear all;
load SkSim.mat

%################## Synthetic signal to work with #########################
% Design parameters
n0 = 0;         % sample asociated to the initial time 
N = length(SVec);        % Number of accumulations
Ta = tVec(2) - tVec(1);      % Accumulation interval

Z = SVec;
%###################### Implementation of the paper #######################
% -------------------------------------------------------------------------
% Paper: "Single-Tone Parameter Estimation from  Discrete-Time
% Observations" by David C. Rife.
%
%   Zn = Xn + j*Yn
%
% where:
%
%   Xn = s(tn) + W(tn),     0 <= n <= N-1
%   Yn = s_(tn) + W_(tn),   0 <= n <= N-1
%
% s  = b0 * cos(omega0*t + theta0)
% s_ = b0 * sin(omega0*t + theta0)
%
% W  is independent gaussian noise with zero mean and variance sigma^2
% W_ is independent gaussian noise with zero mean and variance sigma^2
%
% We assume the first sample taken at t = t0 = n0*T
%--------------------------------------------------------------------------

% -----------------------------  Bounds  ----------------------------------
sigma = 1;% standard deviation of the noise W

% Fisher information matrix
b0 = 4;     % amplitude of the real/imag sinusoidal
n0 = 0;     % how many sampling periods (T) elapsed until t0
T =  0.01;  % sampling period
N =  100;   % number of samples involved in the estimation 
P = N*(N-1)/2;
Q = N*(N-1)*(2*N-1)/6;

J = [b0^2 * T^2 * (n0^2 * N + 2*n0*P + Q), 0, b0^2 * T * (n0*N + P);...
                    0                    , N,           0;...
             b0^2*T*(n0*N + P)           , 0,          b0^2 * N] / sigma^2;

Jinv = inv(J);

% compute the Cramer-Rao Lower Bound
CRLB_omega = Jinv(1,1)   % ====> 0.0075
CRLB_b = Jinv(2,2)       % ====> 0.01
CRLB_theta = Jinv(3,3)   % ====> 0.0024

% ------------------------  Maximum-likelihood  ---------------------------
% It can be obtained as:
%
%       L = 2*b* Re[exp(-j*theta)*exp(-j*omega*t0)*A(omega)] - b^2 
%, 
%
%       A(omega) = (1/N) * sum( Zn*exp(-j*n*omega*T) )
%
% Let's say all three parameters [omega, b, theta] are unknown, then:
%
% solve for 
%
%   omegaML = argmax |A(omegaML)|^2
% 
% Then,
%
%   bML = |A(omegaML)|
%
% Finally,
%
%   thetaML = imag[exp(-j*omegaML*t0)*A(omegaML)]
%
%--------------------------------------------------------------------------
nfft = 2^(nextpow2(length(Z))+8);
A  = fft(Z, nfft) / N;
Af = (0:nfft-1)/(Ta*N)*N/nfft;

[b_hat, idx_Amax] = max(abs(A));

f_hat       = Af(idx_Amax)           % ====> 29.2
b_hat                                % ====> 11.0
theta_hat   = angle(A(idx_Amax))     % ====> 0.4
