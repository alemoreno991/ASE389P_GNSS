clc; close all; clear all
addpath('preprocessing')
addpath("acquisition")

%--------------------------------------------------------------------------
% --------------------    Load input signal     ---------------------------
% --------------------            &             ---------------------------
% -------------------- represent it in baseb-and ---------------------------
%--------------------------------------------------------------------------
display('Loading data...')
filename = 'rawintegersamples_fe.bin'; % 
fc = 1575.42e6;
Tfull   = 60;          % Time interval of data to load (1 ms = 1 data chip).
Fs      = 9.6e6;     % Sampling frequency (Hz).
bandpass.flag = true; % The file has the signal in bandpass representation
bandpass.fIF = 2.391428571429e6; % Intermediate freq. of the bandpass

input = InputSignal(filename, Tfull, Fs, bandpass, 'low');
[tau, xVec, Ts, sMix] = input.getBasebandRepresentation();

%%
%--------------------------------------------------------------------------
% --------------------         Acquisition      ---------------------------
%--------------------------------------------------------------------------
cfgAcquisition.method = 'FFT';
cfgAcquisition.fDVtr  = -4000:100:4000;
cfgAcquisition.nStep  = 1;
cfgAcquisition.nTc    = 10;
cfgAcquisition.nAccum = 1;

rcv = Acquisition(tau, xVec, Ts, cfgAcquisition);
for TXID = 8:8 %[8, 10, 15, 18, 20, 23, 24, 27, 32]
    display(['Acquiring PRN ', num2str(TXID)])
    result{TXID} = rcv.acquireFine(tau, xVec, Ts, TXID, cfgAcquisition);
end

sigmaIQ = result{TXID}.sigmaIQ; % This is representative of the noise. (not dependent of the TXID)
varIQ = sigmaIQ^2;

% plotSk(TXID, tau, result, cfgAcquisition); % DEBUGGING
%%
%--------------------------------------------------------------------------
% --------------------          Tracking       ----------------------------
%--------------------------------------------------------------------------
% Get the bandpass representation of the input signal
[tau, xVec, Ts, sMix] = input.getBandpassRepresentation();

% Configure the tracker
cfgTrack.fc     = fc;
cfgTrack.fIF    = bandpass.fIF;
cfgTrack.Ts     = Ts;
cfgTrack.Nc     = 1023;
cfgTrack.Tc     = 1e-3/cfgTrack.Nc;
cfgTrack.nTc    = 1;
cfgTrack.Ta     = cfgTrack.nTc * cfgTrack.Tc *cfgTrack.Nc;
cfgTrack.sMix   = sMix;

cfgTrack.bufferSk.lenght = 100;
cfgTrack.sigmaIQ         = sigmaIQ;

cfgTrack.pll.Bn     = 10;
cfgTrack.pll.order  = 3;

cfgTrack.dll.Bn     = 0.1;

SV = {};
for TXID = 8:8 %[8, 10, 15, 18, 20, 23, 24, 27, 32]
    display(['Tracking PRN ', num2str(TXID)])
    tracker = Tracker(result{TXID}, TXID, cfgTrack);
    % Run the feedback tracking loop
    theta_hat = [];
    fDk_hat = [];
    tsk_hat = [];
    SkdB = [];
    Sk = [];
    CN0_hat = [];
    t = [];

    % Align the 0-th accumulation interval to the code to avoid data bit
    % transition problems (At least when using 1ms intervals)
    Nk = cfgTrack.Ta/cfgTrack.Ts;
    jk0 = floor(mod(result{TXID}.tsk_hat, 1e-3)/cfgTrack.Ts);
    while jk0 < length(xVec) - Nk
        jVeck = (floor(jk0):floor(jk0+Nk-1))';
        tVeck = jVeck*Ts;
        xVeck = xVec(jVeck+1); % Signal for the k-th accumulation interval
        track_result = tracker.update(tVeck, xVeck);
    
        % Save stuff for plotting
        fDk_hat = [fDk_hat track_result.fD_hat];
        theta_hat = [theta_hat track_result.theta_hat];
        tsk_hat = [tsk_hat track_result.tsk_hat + result{TXID}.tsk_hat];
        SkdB = [SkdB track_result.SkdB];
        Sk = [Sk track_result.Sk];
        CN0_hat = [CN0_hat track_result.CN0];
        t = [t tau(floor(jk0+Nk-1))];
    
        % Prepare for next index iteration
        jk0 = jk0 + Nk;
    end

    SV{TXID}.fDk_hat = fDk_hat;
    SV{TXID}.theta_hat = theta_hat;
    SV{TXID}.tsk_hat = tsk_hat;
    SV{TXID}.SkdB = SkdB;
    SV{TXID}.Sk = Sk;
    SV{TXID}.CN0_hat = CN0_hat;
    SV{TXID}.t = t;
end

display('READY')
%% --------------------------- Plots ------------------------------------%%
clc; close all; clear all
load("Tracking.mat");


for TXID = [8, 10, 15, 18, 20, 23, 24, 27, 32]
    analizeResults(TXID, SV);
end