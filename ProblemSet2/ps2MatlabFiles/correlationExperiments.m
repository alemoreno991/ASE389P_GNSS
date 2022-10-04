% correlationExperiments.m
%
% Experiment with properties of pseudorandom sequences.


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
codeType = 'pi';

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

%----- Compute the sequence autocorrelation
[Rseq1,iiVecSeq] = ccorr(X1,X1);
[Rseq2,iiVecSeq] = ccorr(X2,X2);

%----- Compute the sequence crosscorrelation
[Rseq12,iiVecSeq] = ccorr(X1,X2);

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
