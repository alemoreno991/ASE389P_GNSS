clc; close all; clear all;
addpath('preprocessing')
addpath("acquisition")

%--------------------------------------------------------------------------
% --------------------    Load input signal     ---------------------------
% --------------------            &             ---------------------------
% -------------------- represent it in baseband ---------------------------
%--------------------------------------------------------------------------
filename = 'dataout_raw_trimmed_158.bin';
Tfull   = 1;          % Time interval of data to load (1 ms = 1 data chip).
Fs      = 40e6/7;     % Sampling frequency (Hz).
bandpass.flag = true; % The file has the signal in bandpass representation
bandpass.fIF = 1.405396825396879e6; % Intermediate freq. of the bandpass

input = InputSignal(filename, Tfull, Fs, bandpass, 'high');
[tau, xVec, Ts, sMix] = input.getBasebandRepresentation();

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

%--------------------------------------------------------------------------
% --------------------          Tracking       ----------------------------
%--------------------------------------------------------------------------
%% Initialize.
% |Sk|^2 running average buffer.
buffer_length = 100;
Sk2buffer = max(txid14_sk2(:)) * acqparam.varIQ * ones(1, buffer_length); % TODO: Maybe don't make it so I have to multiply by varIQ...

% Phase loop filter parameters.
Bn_target = 10;
loop_order = 1;
[Ad, Bd, Cd, Dd, Bn_act] = configureLoopFilter(Bn_target, Taccum, loop_order);
plls.Ad = Ad;
plls.Bd = Bd;
plls.Cd = Cd;
plls.Dd = Dd;

% Initial phase loop filter state.
if loop_order > 1
    [v,d] = eig(Ad);
    ind = find(diag(d)==1,1);
    xk = v(:,ind);
    xk = 2*pi*txid14_doppler/dot(Cd,xk) * xk;
    xkp1 = xk;
else
    xk   = zeros(0,1);
    xkp1 = zeros(0,1);
end

% Delay loop parameters.
Sm             = 1; % TODO: Make this change on low/high side mixing.
dlls.Bn_target = Bn_target;
dlls.IsqQsqAvg = mean(Sk2buffer);
dlls.sigmaIQ   = sqrt(acqparam.varIQ);
dlls.vp        = 0;
dlls.Tc        = Tc;
dlls.Ip        = 0;
dlls.Qp        = 0;
dlls.Ie        = 0;
dlls.Qe        = 0;
dlls.Il        = 0;
dlls.Ql        = 0;

% Correlator parameters and data.
corrparam.Fs      = Fs;             % sampling frequency, Hz
corrparam.Fif     = fIF;            % intermediate frequency, Hz
corrparam.eml     = 0.5*Tc;         % delay between early and late taps on code, sec
corrparam.txid    = 14;             % TXID/PRN number
corrparam.Taccum  = Taccum;         % accumulation period
corrparam.Tc      = Tc;             % chipping period (1 ms/1023 chips)
corrparam.init_th = 0;              % initial phase estimate
corrparam.init_ts = txid14_tstart;  % initial code start time estimate
corrdata          = [];             % Note the correlator data will initialize in the data.

% Initial estimate of phase, code start time rates.
vcode_k  = 0;
vtheta_k = 0;

% Trigger for Taccum intervals.
% TODO: Integrate this into the correlator data? Maybe. Maybe not.
elap_accums_last = 0;

% Recorders.
record.phase_est  = nan([1 length(signal)]);
record.tstart_est = nan([1 length(signal)]);
record.vtheta_k   = nan([1 length(signal)]);
record.vcode_k    = nan([1 length(signal)]);
record.Sk2        = nan([1 length(signal)]);
record.Skepl      = nan([3 length(signal)]);

%% Loop.

% Use fake signal
Ts = 1/Fs;
taujv = Ts * (0:(length(signal)-1));

Tmax  = 0.2;
pad_l = round(0.15/Ts);

slope  = 2*pi;
th_tru = [slope*taujv(1:pad_l).*ones(1,pad_l), slope*taujv(pad_l)*ones(1, length(taujv)-pad_l)];

signal = cos(2*pi*fIF*taujv+th_tru); % signal, sqrt(W)
corrparam.init_th = 0;          % initial phase estimate
corrparam.init_ts = 0;          % initial code start time estimate
Sk2buffer = [0 0];

xk   = zeros(size(xk));
xkp1 = xk;

tic
% For every value of the signal...  
% for j = 1:length(signal)
for j = 1:round(Tmax/Ts)

    if mod(j,round(.01*Tmax/Ts))==0
        toc
        tic
        disp(round(j/round(Tmax/Ts)*100))
    end

    xj = signal(j);

    % Pass into correlator.
    [Skepl, corrdata] = correlator(xj, j, vtheta_k, vcode_k, corrdata, corrparam);

    % Detect Taccum intervals. If the correlator's n_accums has
    % incremented, that indicates we have a fresh value of Sk and should
    % update the code/phase rate estimates.
    if corrdata.elap_accums ~= elap_accums_last

        % Update the moving window average of |Sk|^2.
        Sk2buffer(2:end) = Sk2buffer(1:end-1);
        Sk2buffer(1)     = abs(Skepl(2)).^2;

        % Carrier loop filter.
        plls.Ip = real(Skepl(2));
        plls.Qp = imag(Skepl(2));
        plls.xk = xkp1;
        [xkp1, vtheta_k] = updatePll(plls);

        % Code loop filter, applying Sm/wc gain.
        dlls.vp = vtheta_k * Sm/(2*pi*fcL1);
        dlls.IsqQsqAvg = mean(Sk2buffer);

        dlls.Ie = real(Skepl(1)); dlls.Qe = imag(Skepl(1));
        dlls.Ip = real(Skepl(2)); dlls.Qp = imag(Skepl(2));
        dlls.Il = real(Skepl(3)); dlls.Ql = imag(Skepl(3));
        
        v_codek = updateDll(dlls);

    end

    % Remember elapsed accumulations for triggering.
    elap_accums_last = corrdata.elap_accums;

    % Record data.
    record.phase_est (j) = corrdata.phase_est;
    record.tstart_est(j) = corrdata.tstart_est;
    record.vtheta_k  (j) = vtheta_k;
    record.vcode_k   (j) = vcode_k;
    record.Skepl   (:,j) = Skepl;

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