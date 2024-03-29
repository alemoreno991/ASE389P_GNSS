function [fdk_hat, tsk_hat, isAcquired] = acquireGPS(tau, x, T, nTc,...
                                                fStep, nStep, TXID, method)
% Spreading code replica
Nc = 1023;
Tc = 1e-3/Nc;
code = generatePRN(TXID); 

% generate data for the k-th accumulation 
% Nk = floor(nTc*(Nc*Tc)/T);
% Ta = Nc*Tc;
Ta = nTc*Nc*Tc;
Nk = floor(Ta/T);
tau = tau(1:Nk);
x   = x(1:Nk);

% resample C at TsamplingIQ
code = repmat(code, [nTc 1] );
C = oversampleSpreadingCode(code, T/Tc, 0, Nk, Nc);

epsilon = 50*(nTc/Ta);
fDk_hat = -epsilon:fStep:epsilon;
M = zeros(Nk,length(fDk_hat));
if strcmp(method, 'FFT')
    Cr = fft(C);
    parfor fDk_idx = 1:length(fDk_hat)
        % phase estimate over the accumulation interval
        theta_jk_hat = 0;
        theta_hat = 2*pi*fDk_hat(fDk_idx)*tau + theta_jk_hat;

        x_tilde = x.*exp(-i*theta_hat);
        Xr_tilde = fft(x_tilde);
        Zr = Xr_tilde.*conj(Cr);
        zk = ifft(Zr);
        
        M(:,fDk_idx) = abs(zk).^2;
    end
else
    parfor fDk_idx = 1:length(fDk_hat)
        for n = 1:nStep:Nk            
            % phase estimate over the accumulation interval
            theta_jk_hat = 0;
            theta_hat = 2*pi*fDk_hat(fDk_idx)*tau + theta_jk_hat;
            
            % Align C with the code in the incoming data.
            Cshift = circshift(C,n-1);

            % Calculate Sk and store in a matrix that represents the 2D-grid
            Sk = sum(x.*exp(-i*(theta_hat)).*Cshift);
            M(n, fDk_idx) = abs(Sk)^2;
        end
    end
end

% % DEBUG
tau_jk = tau(1:Nk)*1e6;
ii = find(tau_jk<1000);
figure();
h = surf(fDk_hat, tau_jk(ii), M(ii,:));
set(h,'LineStyle','none')
title('2D grid - S_k')
xlabel('$\hat{f}_{Dk} [Hz]$','Interpreter','latex')
ylabel('$\tau_{jk} [us]$','Interpreter','latex')

% Find the maximum Sk in the 2D-grid
[max_Sk, max_idx] = max(M(:));
[n, fDk_idx]=ind2sub(size(M),max_idx);

% Estimates
isAcquired = max_Sk > nTc*3e9;
fdk_hat = fDk_hat(fDk_idx); % Hz
tsk_hat = tau(n)*1e6;      % us 



end

