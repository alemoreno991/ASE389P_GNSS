clc; close all; clear all;
addpath('preprocessing')
addpath("acquisition")

%--------------------------------------------------------------------------
% --------------------    Load input signal     ---------------------------
% --------------------            &             ---------------------------
% -------------------- represent it in baseband ---------------------------
%--------------------------------------------------------------------------
filename = 'dfDataHead.bin';
fc = 1575.42e6;
Tfull   = 1;          % Time interval of data to load (1 ms = 1 data chip).
Fs      = 40e6/7;     % Sampling frequency (Hz).
bandpass.flag = true; % The file has the signal in bandpass representation
bandpass.fIF = 1.405396825396879e6; % Intermediate freq. of the bandpass

input = InputSignal(filename, Tfull, Fs, bandpass, 'high');
[tau, xVec, Ts, sMix] = input.getBasebandRepresentation();

%--------------------------------------------------------------------------
% --------------------         Acquisition      ---------------------------
%--------------------------------------------------------------------------
cfgAcquisition.method = 'FFT';
cfgAcquisition.fDVtr  = -50000:100:0;
cfgAcquisition.nStep  = 1;
cfgAcquisition.nTc    = 10;
cfgAcquisition.nAccum = 2;

detected = [];
fDk_hat = [];
tsk_hat = [];
CN0 = [];
rcv = Acquisition(tau, xVec, Ts, cfgAcquisition);
for TXID = 1:32
    result{TXID} = rcv.acquireFine(tau, xVec, Ts, TXID, cfgAcquisition);
end

%--------------------------------------------------------------------------
% --------------------          Tracking       ----------------------------
%--------------------------------------------------------------------------
cfgTrack.fc     = fc;
cfgTrack.fIF    = bandpass.fIF;
cfgTrack.Fs     = Fs;
cfgTrack.Tc     = 1e-3;
cfgTrack.Ta     = cfgAcquisition.nTc * cfgTrack.Tc;
cfgTrack.sMix   = sMix;

cfgTrack.bufferSk.lenght = 100;
cfgTrack.sigmaIQ         = result{TXID}.sigmaIQ;

cfgTrack.pll.Bn     = 10;
cfgTrack.pll.order  = 3;

cfgTrack.dll.Bn     = 0.01;

tracker = Tracker(estimationSV, TXID, cfg);
% Run the feedback tracking loop
for jk= jk0:Nk:length(tau)
    xVeck = xVec(jk:jk+Nk-1); % Signal for the k-th accumulation interval
    result = tracker.update(xVeck);
end

%% plot bla
% Plot fake signal
figure(2), clf
subplot(211)
hold on, grid on
plot(taujv, signal)
xlim([0 5/fIF])
xlabel('Time (s)')
ylabel('Amplitude')
title('Signal (One Period)')

subplot(212)
hold on, grid on
plot(taujv, th_tru)
plot(taujv, record.phase_est, '--')
xlim([0 Tmax])

figure(3), clf
subplot(121)
hold on, grid on
plot(taujv, record.vtheta_k)

subplot(122)
hold on, grid on
scatter(real(record.Skepl(2,:)), imag(record.Skepl(2,:)))
axis equal