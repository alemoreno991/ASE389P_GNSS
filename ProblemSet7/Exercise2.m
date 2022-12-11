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
Tfull   = 10;          % Time interval of data to load (1 ms = 1 data chip).
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

% Wenkai's load data.
Toffset = 0;                    % Time from beginning to start taking samples.
% Tfull   = Tfull;                % Time interval of data to load (1 ms = 1 data chip).
Fs      = 40e6/7;               % Sampling frequency (Hz).
Ts      = 1/Fs;                 % Sampling period (s).
Noffset = floor(Fs*Toffset);    % Number of data samples to offset start.
N       = floor(Fs*Tfull);      % Number of data samples to load.
N       = N+8-(mod(N,8));       % Round to nearest 8.
nfft    = 2^9;                  % Size of FFT used in power spectrum estimation.
fIF     = 1.405396825396879e6;  % Intermediate frequency.

fid           = fopen('dfDataHead.bin','r','l');
[yhist,count] = binloadSamples(fid,N,'dual');
yhist         = yhist(:,1);
fclose(fid);

xVec = yhist;
tVec = (0:length(yhist))*Ts;

%--------------------------------------------------------------------------
% --------------------          Tracking       ----------------------------
%--------------------------------------------------------------------------
cfgTrack.fc     = fc;
cfgTrack.fIF    = fIF;
cfgTrack.Ts     = Ts;
cfgTrack.Nc     = 1023;
cfgTrack.Tc     = 1e-3/cfgTrack.Nc;
cfgTrack.nTc    = 1;
cfgTrack.Ta     = cfgTrack.nTc * cfgTrack.Tc * cfgTrack.Nc;
cfgTrack.sMix   = sMix;

cfgTrack.bufferSk.lenght = 100;
cfgTrack.sigmaIQ         = result{TXID}.sigmaIQ;

cfgTrack.pll.Bn     = 10;
cfgTrack.pll.order  = 3;

cfgTrack.dll.Bn     = 0.1;

tracker = Tracker(result{TXID}, TXID, cfgTrack);
% Run the feedback tracking loop
jk0 = 1;
Nk = cfgTrack.nTc * cfgTrack.Tc *cfgTrack.Nc/Ts;
theta_hat = [];
fDk_hat = [];
tsk_hat = [];
SkdB = [];
Sk = [];
t = [];
n_accums = Tfull/Ts * 1/Nk;

NkOffset = floor(mod(result{TXID}.tsk_hat,1e-3)/cfgTrack.Ts);
for i_accum = 0:n_accums-2
    Nk0 = floor(i_accum*Nk+NkOffset);
    Nkf = floor((i_accum+1)*Nk+NkOffset);
    jVeck = (Nk0:Nkf)';
    jVeck(jVeck>length(xVec)-1) = [];
    tVeck = jVeck*Ts;
    xVeck = xVec(jVeck+1); % Signal for the k-th accumulation interval
    track_result = tracker.update(tVeck, xVeck);

    % Save stuff for plotting
    fDk_hat = [fDk_hat track_result.fD_hat]; %#ok<*AGROW> 
    theta_hat = [theta_hat track_result.theta_hat];
    tsk_hat = [tsk_hat track_result.tsk_hat];
    SkdB = [SkdB track_result.SkdB];
    Sk = [Sk track_result.Sk];
    t = [t tVec(Nkf)];   
end

%% --------------------------- Plots ------------------------------------%%
close all

%------ |Sk|^2 
figure(1), clf
plot(t, SkdB, LineWidth=2)
title('$|Sk|^2$', Interpreter='latex')
xlabel('time [s]')
ylabel('$|Sk|^2$ [dB]', Interpreter='latex')

figure(2), clf
scatter3(real(Sk), imag(Sk), t)
title('Sk', Interpreter='latex')
xlabel('Real')
ylabel('Imag', Interpreter='latex')
lims = [-3000 3000];
xlim(lims)
ylim(lims)

%------ Doppler 
figure(3), clf
plot(t, round(fDk_hat,2), LineWidth=2)
title('Doppler')
xlabel('time [s]')
ylabel('$\hat{f}_D$ [Hz]', Interpreter='latex')

%------ Phase delay
figure(4), clf
plot(t, unwrap(theta_hat), LineWidth=2)
title('Phase delay')
xlabel('time [s]')
ylabel('$\hat{\theta}$ [deg]', Interpreter='latex')

%------ Code delay 
figure(5), clf
plot(t, round(tsk_hat*1e6,2), LineWidth=2)
title('Code delay')
xlabel('time [s]')
ylabel('$\hat{t}_{sk}$ [us]', Interpreter='latex')