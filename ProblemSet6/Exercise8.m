clc; close all; clear all

% Configure the time vector
Ts =   1.7e-6; % sampling rate [sec]
tf = 3;        % final time [sec]

% Configure the emulation of the incoming signal.
Ta     =    10e-3; % accumulation time, sec
Ampl   =       40; % signal amplitude, sqrt(W)
fIF    =    2.5e3; % intermediate frequency, Hz
CN0    =       50; % carrier-to-noise ratio, dB-Hz
theta0   =      0; % phase offset initial value

% Configure the loop filter
Bn    = 10; % target bandwidth for close-loop
order =  2; % order of loop filter D[z]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% Don't modify beyond this point unless you wanna suffer %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instantiate the time vector
tauj = 0:Ts:tf; 
N  = length(tauj); % Total number of samples in the signal

% Number of samples in accumulation interval
Nk = Ta/Ts;        

% Instantiate the incoming signal
[x_tauj, theta] = simulateSignal(tauj, Ampl, CN0, fIF, Ta, theta0);

% Instantiate the loop filter
vk         =                  0; % estimated phase rate [rad]
xkp1       =   zeros(order-1,1); % loop filter state.
[s.Ad,s.Bd,s.Cd,s.Dd,Bn_act] = configureLoopFilter(Bn, Ta, order);

% Auxiliary stuff to implement the feedback control loop for carrier phase
% tracking
Sk         = 0;          % signal accumulator state
theta_hat  = zeros(1,N); % estimated phase time history [rad]

%---------------------- Run the feedback loop -----------------------------
SkVtr = [];
for j = 1:N-1
    if isAccumulationDone(tauj(j), Ta)
        % single update step of a phase tracking loop
        s.Ip = real(Sk);
        s.Qp = imag(Sk);
        s.xk = xkp1;
        [xkp1,vk] = updatePll(s);
        
        % Save Sk for post-processing
        SkVtr = [SkVtr, Sk];

        % Dump integral.
        Sk = 0;
    end

    % Generate local replica.
    r_tauj = exp(-1j*(2*pi*fIF*tauj(j) + theta_hat(j)));

    % Increment accumulator.
    Sk = Sk + x_tauj(j)*r_tauj;

    % Save the currently estimated phase offset.
    theta_hat(j+1) = theta_hat(j) + vk*Ts;
end

%% Plot.
figure(1), clf
subplot(311)
hold on, grid on
plot(tauj, x_tauj, 'k', LineWidth=2)
xlim([0 1/fIF])
xlabel('Time (s)')
ylabel('Amplitude')
title('Signal')

subplot(312)
hold on, grid on
plot(tauj, theta, 'k', LineWidth=2)
plot(tauj, theta_hat, 'r--', LineWidth=2)
xlabel('Time (s)')
ylabel('Phase (rad)')
legend('True','Estimated')
title("Phase Tracking (loop filter order = " + num2str(order) + ")")

subplot(313)
grid on
plot(tauj, theta_hat-theta, 'k', LineWidth=2)
xlabel('Time (s)')
ylabel('Phase (rad)')
legend('Error')
title("Phase Tracking Error (loop filter order = " + num2str(order) + ")")

ESk2 = mean(abs(SkVtr).^2);
varIQ_exp = ESk2/(2*(10^(CN0/10)*Ta+1));