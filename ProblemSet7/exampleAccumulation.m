% This script is a minimalist realization of a single correlation and
% accumulation operation performed on one GPS L1 C/A code.  xVec holds samples
% from the data set dfDataHead.bin.

% PRN for target satellite
txId = 31;
% Intermediate frequency in Hz
fIF = 1.405396825397e6;
% Sampling time interval in seconds
T = 7/40e6; 
% Nominal chip interval in seconds
Tc = 0.001/1023;
% Approximate Doppler (taken from GRID output for PRN 31)
fD = -256.52;
% The Doppler that acquisition and tracking see is opposite fD due to
% high-side mixing
fD_internal = -fD;
% Nk is the number of samples is one 1-ms accumulation.  It's ok for this
% number to be approximate
Nk = 5714;
% One of the code start times is very near the following sample index
jk = 113593;
% Time vector covering the accumulation
tVec = [0:Nk-1]'*T;
% Generate the +/-1-valued code (not yet oversampled)
cacode = 2*cacodegn(txId) - 1;     
% Oversample the code
cacode_oversampled_prompt = oversampleSpreadingCode(cacode,T/Tc,0,Nk,1023);
% Generate the phase argument of the local carrier replica
ThetaVec = [2*pi*(fIF + fD_internal)*tVec];
% Generate the local carrier replica
carrierVec = exp(-i*ThetaVec);
% Generate the full local replica, with both code and carrier
lVeck = carrierVec.*cacode_oversampled_prompt;
% Isolate the kth code interval from the data. xVec here holds the +/-1 and
% +/-3-valued data samples from dfDataHead.bin.  The first element in xVec
% holds the first sample in dfDataHead.bin.
xVeck = xVec(jk:jk+Nk-1);
% Perform correlation and accumulation
Sk = sum(xVeck.*lVeck)
% Examine the squared magnitude of Sk in dB.  This should be close to 68.29
% dB 
SkdB = 10*log10(abs(Sk)^2)