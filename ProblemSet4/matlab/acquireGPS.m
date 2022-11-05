function [fdk_hat, tsk_hat, isAcquired] = acquireGPS(tau, x, T, ...
                                                fStep, nStep, TXID, method)
% Spreading code replica
Nc = 1023;
Tc = 1e-3/Nc;
code = generatePRN(TXID); 

% generate data for the k-th accumulation 
Nk = floor(1e-3/T);
Ta = Nk*T;
tau = tau(1:Nk);
x   = x(1:Nk);

% resample C at TsamplingIQ
C = oversampleSpreadingCode(code, Nc/Nk, 0, Nk, Nc);

epsilon = 5*(1/Ta);
fDk_hat = -epsilon:fStep:epsilon;
M = zeros(Nk,length(fDk_hat));
if strcmp(method, 'FFT')
    Cr = fft(C);
    
    % k-th accumulation interval
    j = n + (0:1:Nk-1);
    parfor fDk_idx = 1:length(fDk_hat)
            x_tilde = x(j).*exp(-i*2*pi*fDk_hat(fDk_idx)*tau(j));
            Xr_tilde = fft(x_tilde);
            Zr = Xr_tilde*Cr';
            zk = ifft(Zr);
            
            M(:,fDk_idx) = abs(zk)^2;
%             if max_Sk > 3e9
%                 isAcquired = 1;
%                 fdk_hat = fDk_hat(fDk_idx); % Hz
%                 tsk_hat = tau(k_max)*1e6;   % us 
%             end
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

% DEBUG
tau_jk = tau(1:Nk)*1e6;
h = surf(fDk_hat, tau_jk, M);
set(h,'LineStyle','none')
title('2D grid - S_k')
xlabel('$\hat{f}_{Dk} [Hz]$','Interpreter','latex')
ylabel('$\tau_{jk} [us]$','Interpreter','latex')

% Find the maximum Sk in the 2D-grid
[max_Sk, max_idx] = max(M(:));
[n, fDk_idx]=ind2sub(size(M),max_idx);

% Estimates
isAcquired = max_Sk > 3e9;
fdk_hat = fDk_hat(fDk_idx); % Hz
tsk_hat = tau(n)*1e6;      % us 

end

