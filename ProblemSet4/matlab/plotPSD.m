function [] = plotPSD(Syy, fVec, nfft, fs, yLow, yHigh)

%----- Plot results
T = nfft/fs;
delf = 1/T;
fcenter = (nfft/2)*delf;
fVec = fVec - fcenter;
Syy = [Syy(nfft/2 + 1 : end); Syy(1:nfft/2)];
area(fVec/1e6,10*log10(Syy),yLow);
ylim([yLow,yHigh]);
grid on;
shg;
xlabel('Frequency (MHz)');
ylabel('Power density (dB/Hz)');
%figset
title('Power spectral density estimate of GPS L1 Signal');
shg;

end

