% gridPSD.m
%
% Compute power spectrum of GPS samples taken off the GRID receiver

clear; clc;
addpath("research/toolbox/")
%----- Setup
Tfull = 0.5;         % Time interval of data to load
fs = 40e6/7;         % Sampling frequency (Hz)
N = fs*Tfull;        
N = floor(N/16)*16;  % Number of data samples to load
nfft = 2^9;          % Size of FFT used in power spectrum estimation

%----- Load data
fid = fopen('dfDataHead.bin', 'r','l'); 
[Y,count] = binloadSamples(fid,N,'dual');
Y = Y(:,1);
if(count ~= N)
  error('Insufficient data');
end

%----- Compute power spectrum estimate
[Syy,fVec] = pwelch(Y,hann(nfft),[],nfft,fs);


%----- Plot results
yLow = -63;
yHigh = -55;
area(fVec/1e6,10*log10(Syy),yLow);
ylim([yLow,yHigh]);
grid on;
shg;
xlabel('Frequency (MHz)');
ylabel('Power density (dB/Hz)');
figset
title('Power spectral density estimate of GPS L1 Signal');
shg;



