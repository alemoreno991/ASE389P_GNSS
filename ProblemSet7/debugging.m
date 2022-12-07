%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            Problem Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all; format long g
addpath('preprocessing')
addpath("acquisition")

% Solution to the Problem
fDk = [...
    1121.6 ...
   -1551.17 ...
   -1222.22 ...
   -2711.30 ...
   -2652.97 ...
    3103.75 ...
    2049.57 ...
    2186.67 ...
    2459.39 ...
   -1424.31 ...
];

rho = [ ...
    20985472.06 ...
    19972502.03 ...
    18551086.07 ...
    19765301.28 ...
    20277962.34 ...
    21700953.46 ...
    23269270.58 ...
    22529297.11 ...
    19915944.89 ...
    19728164.68 ...
];

C_NO = [
  44.5 ...
  47.7 ...
  45.8 ...
  47.4 ...
  45.2 ...
  42.8 ...
  40.6 ...
  41.1 ...
  47.0 ...
  45.8 ...
];

SV_present = [1 2 5 10 12 15 21 24 29 30];

%--------------------------------------------------------------------------
% --------------------    Load input signal     ---------------------------
% --------------------            &             ---------------------------
% -------------------- represent it in baseband ---------------------------
%--------------------------------------------------------------------------
filename = 'niData01head_5MHz.bin';
Tfull   = 0.5;          % Time interval of data to load (1 ms = 1 data chip).
Fs      = 5.0e6;     % Sampling frequency (Hz).
bandpass.flag = false; % The file has the signal in baseband representation

input = InputSignal(filename, Tfull, Fs, bandpass, 'high');
[tau, xVec, Ts, ~] = input.getBasebandRepresentation();

%--------------------------------------------------------------------------
% --------------------         Acquisition      ---------------------------
%--------------------------------------------------------------------------
fStep = 100;
config.method = 'FFT';
config.fDVtr  = -5000:fStep:5000;
config.nStep  = 1;
config.nTc    = 10;
config.nAccum = 2;

detected = [];
fDk_hat = [];
tsk_hat = [];
CN0 = [];
rcv = Acquisition(tau, xVec, Ts, config);
for TXID = 1:32
    result{TXID} = rcv.acquireFine(tau, xVec, Ts, TXID, config);
    detected = [detected result{TXID}.isAcquired];
    fDk_hat = [fDk_hat result{TXID}.fDk_hat];
    tsk_hat = [tsk_hat result{TXID}.tsk_hat];
    CN0 = [CN0 result{TXID}.CN0];
end

if SV_present == find(detected == 1)
    if abs(fDk_hat(SV_present) - fDk) < 1.5*fStep
        display('Succesful detection')
    else
        display('Detection worked but the doppler frequencies are off')
    end

    tsk_h = tsk_hat(SV_present);
    for i = 1:length(rho)-1
        tsk(i) = mod((rho(1) - rho(i+1))/physconst('LightSpeed'), 1e-3);
        tsk_check(i) = mod(tsk_h(1) - tsk_h(i+1), 1e-3); 
    end
    
    fDk_hat(SV_present)'
    [round(tsk',6) round(tsk_check',6)]*1e6
else
    display('Detection failed!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                  Exam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all; format long g
addpath('preprocessing')
addpath("acquisition")

SV_present = [1, 11, 14, 18, 23, 31]

%--------------------------------------------------------------------------
% --------------------    Load input signal     ---------------------------
% --------------------            &             ---------------------------
% -------------------- represent it in baseband ---------------------------
%--------------------------------------------------------------------------
filename = 'dataout_raw_trimmed_158.bin';
Tfull   = 1;          % Time interval of data to load (1 ms = 1 data chip).
Fs      = 40e6/7;     % Sampling frequency (Hz).
bandpass.flag = true; % The file has the signal in bandpass representation
bandpass.fIF = 1610476.187; % Intermediate freq. of the bandpass

input = InputSignal(filename, Tfull, Fs, bandpass, 'high');
[tau, xVec, Ts, ~] = input.getBasebandRepresentation();

%--------------------------------------------------------------------------
% --------------------         Acquisition      ---------------------------
%--------------------------------------------------------------------------
config.method = 'FFT';
config.fDVtr  = -50000:100:0;
config.nStep  = 1;
config.nTc    = 10;
config.nAccum = 2;

detected = [];
fDk_hat = [];
tsk_hat = [];
CN0 = [];
rcv = Acquisition(tau, xVec, Ts, config);
for TXID = 1:32
    result{TXID} = rcv.acquireFine(tau, xVec, Ts, TXID, config);
    detected = [detected result{TXID}.isAcquired];
    fDk_hat = [fDk_hat result{TXID}.fDk_hat];
    tsk_hat = [tsk_hat result{TXID}.tsk_hat];
    CN0 = [CN0 result{TXID}.CN0];
end

SV_detected = find(detected==1)
length(SV_detected)
fDk_hat(SV_detected)
round(mod(tsk_hat(SV_detected)*1e6,1000))
round(CN0(SV_detected),1)


%% DEBUG
TXID = 23;
Sk = result{TXID}.Sk;
figure();
h = surf(config.fDVtr, tau(1:size(Sk,1))*1e6, Sk(:,:));
set(h,'LineStyle','none')
title('2D grid - S_k')
xlabel('$\hat{f}_{Dk} [Hz]$','Interpreter','latex')
ylabel('$\tau_{jk} [us]$','Interpreter','latex')