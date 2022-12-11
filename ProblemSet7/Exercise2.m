clc; close all; clear all;
addpath('preprocessing')
addpath("acquisition")

%--------------------------------------------------------------------------
% --------------------    Load input signal     ---------------------------
% --------------------            &             ---------------------------
% -------------------- represent it in baseb-and ---------------------------
%--------------------------------------------------------------------------
filename = 'dfDataHead.bin';
fc = 1575.42e6;
Tfull   = 5;          % Time interval of data to load (1 ms = 1 data chip).
Fs      = 40e6/7;     % Sampling frequency (Hz).
bandpass.flag = true; % The file has the signal in bandpass representation
bandpass.fIF = 1.405396825396879e6; % Intermediate freq. of the bandpass

input = InputSignal(filename, Tfull, Fs, bandpass, 'high');
[tau, xVec, Ts, sMix] = input.getBasebandRepresentation();

%--------------------------------------------------------------------------
% --------------------         Acquisition      ---------------------------
%--------------------------------------------------------------------------
cfgAcquisition.method = 'FFT';
cfgAcquisition.fDVtr  = -5000:100:5000;
cfgAcquisition.nStep  = 1;
cfgAcquisition.nTc    = 10;
cfgAcquisition.nAccum = 2;

rcv = Acquisition(tau, xVec, Ts, cfgAcquisition);
for TXID = 14:14
    result{TXID} = rcv.acquireFine(tau, xVec, Ts, TXID, cfgAcquisition);
end


%%%%%-------------
cfgAcquisition.nTc = 1;% TODO: fix this later
[tau, xVec, Ts, sMix] = input.getBandpassRepresentation();
%--------------------------------------------------------------------------
% --------------------          Tracking       ----------------------------
%--------------------------------------------------------------------------
cfgTrack.fc     = fc;
cfgTrack.fIF    = bandpass.fIF;
cfgTrack.Ts     = Ts;
cfgTrack.Nc     = 1023;
cfgTrack.Tc     = 1e-3/cfgTrack.Nc;
cfgTrack.nTc    = cfgAcquisition.nTc;
cfgTrack.Ta     = cfgAcquisition.nTc * cfgTrack.Tc * cfgTrack.Nc;
cfgTrack.sMix   = sMix;

cfgTrack.bufferSk.lenght = 100;
cfgTrack.sigmaIQ         = result{TXID}.sigmaIQ;

cfgTrack.pll.Bn     = 10;
cfgTrack.pll.order  = 3;

cfgTrack.dll.Bn     = 0.1;

tracker = Tracker(result{TXID}, TXID, cfgTrack);
% Run the feedback tracking loop
jk0 = floor(cfgAcquisition.nAccum * cfgAcquisition.nTc * cfgTrack.Tc * cfgTrack.Nc/ Ts);
Nk = floor( cfgAcquisition.nTc * cfgTrack.Tc *cfgTrack.Nc / Ts);
theta_hat = [];
fDk_hat = [];
tsk_hat = [];
SkdB = [];
Sk = [];
for jk= jk0:Nk:length(xVec)-Nk
    xVeck = xVec(jk:jk+Nk-1); % Signal for the k-th accumulation interval
    track_result = tracker.update(xVeck);
    fDk_hat = [fDk_hat track_result.fD_hat];
    theta_hat = [theta_hat track_result.theta_hat];
    tsk_hat = [tsk_hat track_result.tsk_hat];
    SkdB = [SkdB track_result.SkdB];
    Sk = [Sk track_result.Sk];
end

%% --------------------------- Plots ------------------------------------%%
close all
t = tau(3*Nk:Nk:end);

%------ |Sk|^2 
figure(), clf
plot(t, SkdB, LineWidth=2)
title('$|Sk|^2$', Interpreter='latex')
xlabel('time [s]')
ylabel('$|Sk|^2$ [dB]', Interpreter='latex')

figure(), clf
scatter3(real(Sk), imag(Sk), t)
title('Sk', Interpreter='latex')
xlabel('Real')
ylabel('Imag', Interpreter='latex')


%------ Doppler 
figure(), clf
plot(t, round(fDk_hat,2), LineWidth=2)
title('Doppler')
xlabel('time [s]')
ylabel('$\hat{f}_D$ [Hz]', Interpreter='latex')

%------ Phase delay
figure(), clf
plot(t, round(rad2deg(wrapToPi(theta_hat)),2), LineWidth=2)
title('Phase delay')
xlabel('time [s]')
ylabel('$\hat{\theta}$ [deg]', Interpreter='latex')

%------ Code delay 
figure(), clf
plot(t, round(tsk_hat*1e6,2), LineWidth=2)
title('Code delay')
xlabel('time [s]')
ylabel('$\hat{t}_{sk}$ [us]', Interpreter='latex')