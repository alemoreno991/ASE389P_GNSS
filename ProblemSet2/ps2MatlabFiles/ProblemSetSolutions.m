%% Problem 1
% The power spectral density of a random binary sequence with chip interval
% Tc is given by 
%                           Sx(f) = Tc*(sinc(f*Tc))^2
%
% From lecture we know that over 90% of the power in Sx(f) lies in the 
% interval -1/Tc < f < 1/Tc
%--------------------------------------------------------------------------
clc; close all; clear all
format short
% What percentage of power lies in -2/Tc < f < 2/Tc ?
Tc = 1;
f = linspace(-2/Tc, 2/Tc, 1e5);
Sx = Tc*(sinc(f*Tc)).^2;
plot(f,Sx);
xlabel('freq [Hz]'); ylabel('Sx'); title('Power Spectral Density')

Sx = @(f) Tc*(sinc(f*Tc)).^2;
power = integral(Sx,-2/Tc, 2/Tc);

% Integral of Sx(f) over -Inf < f < Inf = 1. Then, the percentage is 
total_power = 1;
round(power/total_power, 3)


%% Problem 2
% GPS L1 C/A code has a Tc ~ 1us 
% ==> null-to-null bandwidth of the main lobe of the power spectrum (B1)
% B1 = 2/Tc ~ 2 MHz
%
% Note that the spectrum of Sx within B1 is not filled uniformly with power.
%--------------------------------------------------------------------------
clc; close all; clear all
format short
% Since in class we use a different convention for the fourier transform 
% than the one that MATLAB uses by default. Then, I change it right away
oldVal = sympref('FourierParameters',[1 -(2*sym(pi))]);

Tc = 1e-6; W = 2e6;
syms f f1 tau tau1
Sx = Tc*(sinc(f1*Tc))^2;
Sx_efficient = rectangularPulse(f/W)/W; 

figure()
fplot(Sx,[-2/Tc, 2/Tc], Color='k')
hold on
fplot(Sx_efficient,[-2/Tc, 2/Tc], Color='b')
title('Power spectral density')
xlabel('f')
legend('Sx','Sx_{efficient}')
grid on

% a) Find the autocorrelation function Rx(tau) of the spreading waveform by
% taking the inverse transform of Sx(f)
Rx_efficient = ifourier(Sx_efficient,f,tau);

figure()
fplot(Rx_efficient, [-2*Tc, 2*Tc], Color='b')
title('Autocorrelation of Rx_{efficient}')
xlabel('tau')
grid on

% b) How wide is the main peak in the autocorrelation function Rx (from
% first left to first right zero-crossing)?
right_zero_crossing = vpa(solve(Rx_efficient));
main_peak_width = 2*right_zero_crossing % since the function is symmetric

% c) Compare this width to the width of the peak of the autocorrelation
% function Rx_bar for a random binary sequence with a null-to-null
% bandwidth B1 = W
Rx_bar = triangularPulse(-Tc, Tc, tau1); % from class notes

figure()
fplot(Rx_bar,[-2*Tc, 2*Tc], Color='k')
hold on
fplot(Rx_efficient, [-2*Tc, 2*Tc], Color='b')
title('Autocorrelation') 
xlabel('tau')
legend('Rx','Rx_{efficient}')
grid on

main_peak_width_bar = 2*Tc

% d) The following ratio shows that the width of the main peak of the 
% autocorrelation function of the more efficient waveform (in terms of 
% noise) is more than 3 times wider than the original random binary
% sequence. This means that one would loose precision. 
% 
% My logic is the following, one would want the autocorrelation to be a 
% dirac delta. Then, the achieved precision would be excellent. As far as
% one gets from this ideal situation the worse the precision gets.
ratio = main_peak_width / main_peak_width_bar 

% What about the high-C/N0 scenario vs the low-C/N0 scenario?
% The wider your signal is in the time-domain, the narrower it is in the
% frequency domain. The narrower the frequency domain, the further away
% from a Dirac delta the autocorrelation gets. 
% However, the wider the signal in the time domain the better would be the
% receiver's ability to understand the signal under low C/N0 situations.
% Therefore, a tradeoff exists and there's no right answer. I would say
% that in a low-C/N0 scenario one should try to slow the sequence down as
% much as possible to detect it because detecting it and getting low
% accuracy in the range measurement is better than not detecting at all. In
% a high-C/N0 I would by greedy with respect to accuracy and make the
% spreading code as fast as possible.

% e) The pseudo-random binary sequence was chosen because of its built-in
% ranging capabilities, multiple access characteristic and jamming
% rejection. The reason why it needs to be pseudo-random is because without
% knowing the actual sequence there's no way to actually "lock" to the
% spreading code (no way to have a local replica).

%% Problem 3
% Suppose you buy a fancy low-noise amplifier and crygenic cooling system
% for your GNSS receiver so that your receiver temperature drops to 
%                           TR = 50K
% Assuming an antenna temperature of 100K, calculate the corresponding
% noise density N0 in dBW/Hz
clc, close all, clear all

k = 1.380649e-23; % Boltzmann's constant
TA = 100;       % Antenna temperature
TR = 50;        % Receiver temperature
Ts = TA + TR;   % System temperature
N0 = k*Ts;      % Effective noise density for the system

% Effective noise density for the system expressed in dBW/Hz
N0_dBW_Hz = 10*log10(N0) % Effective noise density [dBW/Hz]

% Now consider how this reduction in noise density gets offset by the
% interference effects of other GNSS signals, due to "multiple-access". 
% Supposse that the desired signal is one of M = 10 received signals, each 
% with received power
%                            PI = -155 dbW
%                            Tc = 1/1023 us
%
% Notes from class: multiple-access interference is characterized as follow
%
%                          I0 = (2/3)*PI*Tc*(M-1)
M = 10;
PI_dBW = -155;
PI = 10^(PI_dBW/10); 
Tc = (1/1023) * 1e-3;
I0 = (2/3)*PI*Tc*(M-1);

I0_dBW_Hz = 10*log10(I0)    % Multiple-access interference [dBW/Hz]

% What will be the equivalent (noise + interference) density N0_eq?
N0eq = N0 + I0;
N0eq_dBW_Hz = 10*log10(N0eq) % Noise equivalent density [dBW/Hz]

%% Problem 4
% Suppose that C(t) has a power spectral density given by:
%
%                   Sc(f) = Tc * rectangularPulse(f*Tc)
%
% Where 'rectangularPulse' is the rect function introduced in lecture.
% Imagine a constellation of GNSS satellites broadcasting signals having
% this new Sc(f). Assume the signals have no data bit modulation, so that
% the signal's baseband representation is of the form:
%
%                   r(t) = sqrt(P)*C(t)*exp[j theta(t)]
%
clc, close all, clear all

% For this case, calculate I0 = SI(f=0) for a single multiple-access
% interference signal with power PI, where I0 and SI(f), the power spectrum
% of I(t), were defined in class.
%
% ANSWER:
% Following the same reasoning as in lecture:
% If CI(t) is a random-binary sequence "like C(t)" (same Tc), and if the
% two are uncorrelated. Then, (assuming fD=0)
%                       
%                   SrI(f) = PI*Sc(f)
% Thus, 
%                   SI(f) = Sc(f) x SrI(f) = Sc(f) x (PI*Sc(f))
%
%                   SI(f) = PI * integral_{-inf}^{inf}(Sc(f-x)Sc(x) dx
%
% Since, Sc(f) = Sc(-f)
%
%                   SI(f) = PI * integral_{-inf}^{inf}(Sc(x-f)Sc(x) dx
%
% Now, using the fact that the receiver has a low pass filter, only the f=0
% part of the spectral density leaks through. Therefore, it's ok to assume 
% I0 = SI(f=0) (REMEMBER: I0 is also a power density)
%
%     SI(0) = PI * integral_{-inf}^{inf}(Sc(x)^2) dx
%
%     SI(0) = PI * integral_{-inf}^{inf}(Tc*rectangularPulse(x*Tc)^2) dx
%
%     SI(0) = PI*Tc
I0 = PI*Tc;
I0_dBW_Hz = 10*log10(I0)

% Now consider M multiple-access signals (including the desired signal). If
% each of the M received signals has power PI, at what value of M does I0
% exceed N0? Express your answer in terms of N0, Tc, PI.
%
% ANSWER:
% The multiple-access interference is characterized by:
%
%           I0 = PI*Tc*(M-1)
%
% We want to know what M would result in I0 >= N0
%
%           PI*Tc*(M-1) >= N0
%           M >= N0/(PI*Tc) + 1, where M is integer
%        
% => I0 exceeds N0 when:
%           M = ceil(N0/(PI*Tc) + 1)
%

%% Problem 5
% TODO: solve this problem (Doesn't look hard. Just plot stuff)

%% Problem 6
clear all; clc; close all;
addpath("research/toolbox/")

%----- Load data
load('prn31_22apr03_01hrs40min00sec_gmt_fl1_46_08mhz_250msec.mat');

start = 8000;
Y_windowed = Y(start:start+400,:);

% Plotting the complex representation of the signal such that it forms a
% "perfect" square in the the complex plane.
figure()
plot(Y_windowed,'.')
xlabel('In phase');
ylabel('Quadrature');
title('Complex representation of the baseband signal');

t = XDelta*(0:400) * 1e6;
CA_code = real(Y_windowed);
Y_code = imag(Y_windowed);

% Plotting the real and imaginary parts of the signal
figure()
subplot(2,1,1)
    plot(t,CA_code)
    xlabel('time [us]');
    ylabel('X_{In-phase}');
    title('Time-domain representation of X_{In-phase}');

subplot(2,1,2)
    plot(t, Y_code)
    xlabel('time [us]');
    ylabel('X_{Quadrature}');
    title('Time-domain representation of X_{Quadrature}');

% TODO: Measure the smallest chip interval for both codes.

% The ratio of the C/A code amplitude to the P(Y) code amplitude
ratio = mean(sqrt(CA_code.^2)) / mean(sqrt(Y_code.^2)) % Amplitude of C/A_code is ~1.4x P(Y) code amplitude 

%% Problem 7
% Verify by "Matlab experiment" that the power spectral density of a random
% binary sequence is given by 
%               
%                   Sx(f) = Tc*sinc(f*Tc)^2
% 
clear all; clc; close all

Tc = 1e-6;
N = 2^16;
X = 2*randi([0 1], N, 1) - 1;
t = (0:N-1)*Tc;

OS = (5/83);
Tos = OS*Tc;
Xos = oversampleSpreadingCode(X, OS, 0, ceil(N/OS), N);
tos = (0:(length(Xos)-1))*Tos;

% Plot the sequence (for better visualization reduce the number of points
figure()
plot( tos, Xos,'o', 'Color', 'red')
hold on
stairs( t, X, '-b', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Signal')
legend('Oversampled', 'Original')

%----- Compute power spectrum estimate
fs = 1/Tos;
nfft = 2^10;
[S_Xos,fVec] = pwelch(Xos,hann(nfft),[],nfft,fs,"twosided");

%----- Plot results
yLow = 10*log10(min(S_Xos));
yHigh = 10*log10(max(S_Xos)) + 10;
T = nfft/fs;
delf = 1/T;
fcenter = (nfft/2)*delf;
fVec = fVec - fcenter;
S_Xos = [S_Xos(nfft/2 + 1 : end); S_Xos(1:nfft/2)];

figure()
area(fVec/1e6, 10*log10(S_Xos), yLow);
ylim([yLow,yHigh]);
grid on;
shg;
xlabel('Frequency (MHz)');
ylabel('Power density (dBW/Hz)');
title('Power spectral density estimate of the oversampled binary random sequence');
shg;

% Height of the main lobe in dBW/Hz 
height_main_lobe = 10*log10(Tc*(sinc(0))^2)

%% Problem 9a

% correlationExperiments.m
%
% Experiment with properties of pseudorandom sequences.


clear;clc;
addpath("research/toolbox/")
%----- Setup
nStages = 10;                 % Number of stages in LFSR
Tc = 1e-3/1023;               % Chip interval in seconds
delChip = 3/217;              % Sampling interval in chips
delOffset  = 0;               % Offset of first sample
delt = delChip*Tc;            % Sampling interval in seconds
fs = 1/delt;                  % Sampling frequency in Hz
Np = 2^nStages - 1;           % Period of the sequence in chips
Nr = 20;                      % Number of repetitions of the sequence
Ns = round(Nr*Np/delChip);    % Number of samples of the sequence 

N = 10000;
varVec = zeros(N,1);
for i=1:N
    % codeType:
    % rand ---- Sequence derived from Matlab randn function
    % pi ------ Sequence derived from the digits of pi
    % mseq ---- Maximal-length sequence with n = nStages
    codeType = 'rand';
    
    %----- Generate codes
    X1 = zeros(Np,1);
    X2 = zeros(Np,1);
    if(strcmp(codeType,'rand'))
      X1 = sign(sign(randn(Np,1)) + 0.1);
      X2 = sign(sign(randn(Np,1)) + 0.1);
    elseif(strcmp(codeType,'pi'))
      [sPi,vPi] = pi2str(2*Np);
      X1 = vPi(1:Np) >= 5;
      X1 = 2*X1 - 1;
      X2 = vPi(Np+1:2*Np) >= 5;
      X2 = 2*X2 - 1;
    elseif(strcmp(codeType,'mseq'))
    ciVec1 = [9, 4]';  
    ciVec2 = [9, 2]';
    a0Vec1 = [1;zeros(nStages-1,1)];
    a0Vec2 = ones(nStages,1);
    X1 = generateLfsrSequence(nStages,ciVec1,a0Vec1);
    X2 = generateLfsrSequence(nStages,ciVec2,a0Vec2);
    X1 = 2*X1 - 1;
    X2 = 2*X2 - 1;
    else
      error('Unrecognized code type');
    end
    
    %----- Compute the sequence crosscorrelation
    [Rseq12,iiVecSeq] = ccorr(X1,X2);
    index = find(iiVecSeq==0);
    crosscorrVec(i) = Rseq12(index);
end

% The problem asks about the variance of the cross-correlation at k=0
variance = var(crosscorrVec)

%% Problem 9b

% correlationExperiments.m
%
% Experiment with properties of pseudorandom sequences.


clear;clc; format short
addpath("research/toolbox/")
%----- Setup
nStages = 10;                 % Number of stages in LFSR
Tc = 1e-3/1023;               % Chip interval in seconds
delChip = 3/217;              % Sampling interval in chips
delOffset  = 0;               % Offset of first sample
delt = delChip*Tc;            % Sampling interval in seconds
fs = 1/delt;                  % Sampling frequency in Hz
Np = 2^nStages - 1;           % Period of the sequence in chips
Nr = 20;                      % Number of repetitions of the sequence
Ns = round(Nr*Np/delChip);    % Number of samples of the sequence 

N = 6;
varVec = zeros(N,1);
for i=1:N
    ciVec1 = [
        [10, 9, 8, 5]',...  
        [10, 9, 7, 6]',...
        [10, 9, 7, 3]',...
        [10, 9, 6, 1]',...
        [10, 9, 5, 2]',...
        [10, 9, 4, 2]'
    ];
    ciVec2 = [
        [10, 9, 8, 7, 5, 4]',...  
        [10, 9, 8, 7, 4, 1]',...
        [10, 9, 8, 7, 3, 2]',...
        [10, 9, 8, 6, 5, 1]',...
        [10, 9, 8, 6, 4, 2]',...
        [10, 9, 8, 6, 4, 2]'  
    ];
    a0Vec1 = [1;zeros(nStages-1,1)];
    a0Vec2 = ones(nStages,1);
    X1 = generateLfsrSequence(nStages,ciVec1(:,i),a0Vec1);
    X2 = generateLfsrSequence(nStages,ciVec2(:,i),a0Vec2);
    X1 = 2*X1 - 1;
    X2 = 2*X2 - 1;
    %----- Compute the sequence crosscorrelation
    [Rseq12,iiVecSeq] = ccorr(X1,X2);
    
    bound_crosscorrVec(i) = max(Rseq12);

    %----- Compute the sequence autocorrelation
    [Rseq11,iiVecSeq] = ccorr(X1,X1);
    index = find(iiVecSeq==0);
    autocorrVec(i) = Rseq11(index);
end

% The value of the autocorr at k=0 ~ (N-1)
autocorrVec

% The bound of the cross-correlation >= sqrt(N)
bound_crosscorrVec

%% Problem 9c

clear;clc;close all
addpath("research/toolbox/")
%----- Setup
nStages = 10;                 % Number of stages in LFSR
Tc = 1e-3/1023;               % Chip interval in seconds
delChip = 3/217;              % Sampling interval in chips
delOffset  = 0;               % Offset of first sample
delt = delChip*Tc;            % Sampling interval in seconds
fs = 1/delt;                  % Sampling frequency in Hz
Np = 2^nStages - 1;           % Period of the sequence in chips
Nr = 20;                      % Number of repetitions of the sequence
Ns = round(Nr*Np/delChip);    % Number of samples of the sequence 
% codeType:
% rand ---- Sequence derived from Matlab randn function
% pi ------ Sequence derived from the digits of pi
% mseq ---- Maximal-length sequence with n = nStages
codeType = 'mseq'; % Change this variable to see the different results

%----- Generate codes
X1 = zeros(Np,1);
X2 = zeros(Np,1);
if(strcmp(codeType,'rand'))
  X1 = sign(sign(randn(Np,1)) + 0.1);
  X2 = sign(sign(randn(Np,1)) + 0.1);
elseif(strcmp(codeType,'pi'))
  [sPi,vPi] = pi2str(2*Np);
  X1 = vPi(1:Np) >= 5;
  X1 = 2*X1 - 1;
  X2 = vPi(Np+1:2*Np) >= 5;
  X2 = 2*X2 - 1;
elseif(strcmp(codeType,'mseq'))
ciVec1 = [9, 4]';  
ciVec2 = [9, 2]';
a0Vec1 = [1;zeros(nStages-1,1)];
a0Vec2 = ones(nStages,1);
X1 = generateLfsrSequence(nStages,ciVec1,a0Vec1);
X2 = generateLfsrSequence(nStages,ciVec2,a0Vec2);
X1 = 2*X1 - 1;
X2 = 2*X2 - 1;
else
  error('Unrecognized code type');
end

%----- Oversample code
X1os = oversampleSpreadingCode(X1,delChip,delOffset,Ns,Np);
X2os = oversampleSpreadingCode(X2,delChip,delOffset,Ns,Np);

%----- Compute autocorrelation 
[R1,iiVec] = ccorr(X1os,X1os);
[R2,iiVec] = ccorr(X2os,X2os);

%----- Compute crosscorrelation 
[R12,iiVec] = ccorr(X1os,X2os);

ratio = max(R1)/max(R12)

%----- Compute power spectra
% The scaling here ensures that sum(Si*delf) = 1 W, as expected, for i = 1, 2.
S1 = abs(delt*fft(X1os)).^2/(Ns*delt);
S2 = abs(delt*fft(X2os)).^2/(Ns*delt);
S12 = abs(delt*fft(R12)).^2/(Ns*delt);
delf = 1/(delt*Ns);
fVec = [0:Ns-1]'*(delf);

%----- Scale for 1-Hz delf for plotting
% We'll present the power spectra in units of dBW/Hz, so we need to scale
% the spectra so that they match a delf = 1 Hz frequency spacing.
S1 = S1*delf;
S2 = S2*delf;
S12 = S12*delf;

%----- Plot
figure(1);clf;
subplot(211)
plot(iiVec,R1/Ns);
grid on;
ylabel('R_{X1}');
title('X1 and X2 autocorrelation')
subplot(212)
plot(iiVec,R2/Ns);
grid on;
ylabel('R_{X2}');
xlabel('Lag (samples)');
figure(2);clf;
plot(iiVec,R12/Ns);
title('X1 and X2 crosscorrelation')
xlabel('Lag (samples)');
grid on;
ylabel('R_{X1,X2}');
figure(3);clf;
subplot(211)
plot(fVec/1e3,10*log10(S1));
grid on;
xlim([0 30]);
ylim([-100,0]);
title('X1 and X2 power spectral densities')
ylabel('S_{X1}(f) (dBW/Hz)');
subplot(212)
plot(fVec/1e3,10*log10(S2));
grid on;
xlim([0 30]);
ylim([-100,0]);
ylabel('S_{X2}(f) (dBW/Hz)');
xlabel('Frequency (kHz)');

%% Extra Problem
clc, close all, clear all
%----- Setup
nStages = 9;                 % Number of stages in LFSR

% f1 (D) =1 + D4 + D9
% f2 (D) =1 + D3 + D5 + D6 + D9
% f3 (D) =1 + D3 + D5 + D8 + D9
% f4 (D) =1 + D3 + D4 + D6 + D9
% f5 (D) =1 + D3 + D7 + D8 + D9
% f6 (D) =1 + D2 + D9
ciVec1 = [9, 4]';  
ciVec2 = [9, 6, 5, 3]';
ciVec3 = [9, 8, 5, 3]';
ciVec4 = [9, 6, 4, 3]';
ciVec5 = [9, 8. 7. 3]';
ciVec6 = [9, 2]';

a0Vec1 = ones(nStages,1);

% Examine the final state to determine if it matches the initial and 
% thus conclude that it's indeed a maximal length LFSR
X1 = generateLfsrSequence(nStages,ciVec1,a0Vec1); % maximal length
X2 = generateLfsrSequence(nStages,ciVec2,a0Vec1); % maximal length
X3 = generateLfsrSequence(nStages,ciVec3,a0Vec1); % not maximal lenght
X4 = generateLfsrSequence(nStages,ciVec4,a0Vec1); % maximal length
X5 = generateLfsrSequence(nStages,ciVec5,a0Vec1); % not maximal length
X6 = generateLfsrSequence(nStages,ciVec6,a0Vec1); % not maximal length

% Now let's see if it's a Gold sequence by analizing the cross-correlation

%----- Compute crosscorrelation 
[R12,iiVec] = ccorr(X1,X2);
[R13,iiVec] = ccorr(X1,X3);
[R23,iiVec] = ccorr(X2,X3);

figure()
subplot(1,3,1)
plot(R12)
title('R12')
grid on

subplot(1,3,2)
plot(R13)
title('R13')
grid on

subplot(1,3,3)
plot(R23)
title('R23')
grid on