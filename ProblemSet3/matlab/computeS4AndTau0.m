function [S4,tau0] = computeS4AndTau0(zkhist,tkhist)
% computeS4AndTau0 : Compute the scintillation index S4 and the decorrelation
%                    time tau0 corresponding to the input complex channel
%                    response function time history zkhist.
%
% INPUTS
%
% zkhist ----- Nt-by-1 vector containing the normalized complex scintillation
%              time history in the form of averages over Ts with sampling
%              interval Ts. zkhist(kp1) is the average over tk to tkp1.
%
% tkhist ----- Nt-by-1 vector of time points corresponding to zkhist.
%
%
% OUTPUTS
%
% S4 --------- Intensity scintillation index of the scintillation time history
%              in zkhist, equal to the mean-normalized standard deviation of
%              the intensity abs(zkhist).^2.
%
% tau0 ------- The decorrelation time of the scintillation time history in
%              zkhist, in seconds.
%
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+
alpha = abs(zkhist);
I = alpha.^2;
S4_squared = (mean(I.^2) - mean(I)^2) / mean(I)^2;
S4 = sqrt(S4_squared);

z_bar = mean(zkhist);
xi = zkhist - z_bar;
[R,lags] = xcorr(xi,'normalized'); 
idx = find(abs(R(find(lags==0):end)) < R(find(lags==0))*exp(-1), 1);
tau0 = tkhist(idx);
