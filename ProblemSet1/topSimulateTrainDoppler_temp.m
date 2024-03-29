% topSimulateTrainDoppler.m
%
% Top-level script for train Doppler simulation


clear; clc; close all;
%----- Setup
fc = 440;
vTrain = 20;
t0 = 0;
x0 = 0;
delt = 0.01;
N = 1000;
vs = 343;
xObs = 56.8;
dObs = 10;

%----- Simulate
[fDVec,tVec] = simulateTrainDoppler(fc,vTrain,t0,x0,xObs,dObs,delt,N,vs);
fApparentVec = fDVec + fc;

%----- Plot
plot(tVec,fDVec + fc, 'r');
xlabel('Time (seconds)');
ylabel('Apparent horn frequency (Hz)');
grid on;
shg;

%----- Generate a sound vector
T = delt*N;                    % simulation time (sec)
fs = 22050;                    % sample frequency (Hz)
deltSamp = 1/fs;               % sampling interval (sec)
Ns = floor(T/deltSamp);        % number of samples
tsamphist = [0:Ns-1]'*deltSamp;
Phihist = zeros(Ns,1);
fApparentVecInterp = interp1(tVec,fApparentVec,tsamphist,'spline');
for ii=2:Ns
  fii = fApparentVecInterp(ii);
  Phihist(ii) = Phihist(ii-1) + 2*pi*fii*deltSamp;
end
soundVec = sin(Phihist);

% %----- Play the sound vector
%sound(soundVec, fs);    

%----- Write to audio file
audiowrite('my_trainout.wav',soundVec,fs);

%----- Write frequency time history to output file
save trainData fApparentVec tVec

%% Extra Points
close all
% I use an spectrogram to better understand the original audio signal
[y, fs] = audioread('trainout.wav');

Nspec = 2^(nextpow2(length(y)) - 7);
wspec = hamming(Nspec);
Noverlap = Nspec/2;
figure()
% spectrogram(y, 512, 128, 128, fs, 'yaxis');
spectrogram(y, wspec, Noverlap, Nspec, fs, 'yaxis');
colormap jet
title("original signal")

set(gca,'Linewidth',1.2,'FontSize',36)
set(gcf,'Position',[2500 100 1550 800])

% Then, I manually analyze my simulated audio signal to make it look as 
% similar to the original as posible by changing the position of the
% receiver
% [y, fs] = audioread('my_trainout.wav');
% 
% figure()
% spectrogram(y, 512, 120, 128, fs, 'yaxis');
% colormap jet
% title("My signal")
% 
% set(gca,'Linewidth',1.2,'FontSize',36)
% set(gcf,'Position',[2500 100 1550 800])