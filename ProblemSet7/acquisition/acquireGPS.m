function [result] = acquireGPS(tau, x, Ts, TXID, config)
% acquireGPS:       This function runs the acquisition algorithm over a
%                   signal `x` to find a specific PRN. In case it succeeds,
%                   the estimated doppler and code-delay are informed back
%                   (among other things).
%
% INPUTS
%
% tau ------------- Time vector 
%
% x   ------------- Signal in complex baseband representation
%
% Ts  ------------- Sampling period 
%
% TXID ------------ ID number of the PRN you want to look for.
%
% config ---------- Structure to configure the acquisition algorithm
% 
%                   
%                   - method: `vanilla` or `FFT` 
%                   
%                   - nStep:   step (number of Ts) used to sweep the 
%                              code-delay axis of the search space (used 
%                              in the vanilla case).
%                   
%                   - fDVtr:   vector representing the doppler axis of
%                              the search space to be swept.
%
%                   - nTc:     number of code periods to be used for the
%                              accumulation interval. (used for each
%                              coherent accumulation)
%
%                   - nAccum:  number of non-coherent accumulations 
%                   
%
% OUTPUTS
%
% info ------------ Structure containing useful information 
%   
%                   - isAcquire: a boolean flag indicating if the TXID was
%                                detected. (ALWAYS check this first. The
%                                values of the rest of the fields in the
%                                structure would be undefined if
%                                `isAcquire = false`)
%
%                   - fdk_hat:   the estimated doppler [Hz]
%
%                   - tsk_hat:   the estimated code-delay [sec]
%
%                   - CN0:       the carrier-to-noise ratio [dB-Hz]
%
%==========================================================================
    % Spreading code replica
    Nc = 1023;
    Tc = 1e-3/Nc;
    code = generatePRN(TXID); 
    
    % generate data for the k-th accumulation 
    Ta = config.nTc*Nc*Tc;
    Nk = floor(Ta/Ts);
    tau = tau(1:Nk);
    x   = x(1:Nk);
    
    % resample C at Ts
    code = repmat(code, [config.nTc 1] );
    C = oversampleSpreadingCode(code, Ts/Tc, 0, Nk, Nc);
    
    % Calculate the Sk matrix for the k-th accumulation
    Sk = calculateSk(tau, x, C, config.fDVtr, config.nStep, Nk, config.method);
     
    % % DEBUG
    % tau_jk = tau(1:Nk)*1e6;
    % ii = find(tau_jk<1000);
    % figure();
    % h = surf(fDk_hat, tau_jk(ii), M(ii,:));
    % set(h,'LineStyle','none')
    % title('2D grid - S_k')
    % xlabel('$\hat{f}_{Dk} [Hz]$','Interpreter','latex')
    % ylabel('$\tau_{jk} [us]$','Interpreter','latex')
    
    % Find the maximum Sk in the 2D-grid
    [max_Sk, max_idx] = max(Sk(:));
    [n, fDk_idx]=ind2sub(size(Sk),max_idx);
    
    % Estimates
    result.isAcquired = max_Sk > 9e9;
    result.fDk_hat = config.fDVtr(fDk_idx); % Hz
    result.tsk_hat = tau(n);           % sec 
    
%     two_sigmaIQ_squared = 61551.7110161536; % mean(M,'all') for satellite 34
    two_sigmaIQ_squared = 568409076.119826; % mean(M,'all') for satellite 34
    result.CN0 = 10*log10((max_Sk - two_sigmaIQ_squared)/(two_sigmaIQ_squared*Ta));

end

function [Sk] = calculateSk(tau, x, C, fDVtr, nStep, Nk, method)
    Sk = zeros(Nk, length(fDVtr));
    if strcmp(method, 'FFT')
        Cr = fft(C);
        for fDk_idx = 1:length(fDVtr)
            % phase estimate over the accumulation interval
            theta_jk_hat = 0;
            theta_hat = 2*pi*fDVtr(fDk_idx)*tau + theta_jk_hat;
    
            x_tilde = x.*exp(-1i*theta_hat);
            Xr_tilde = fft(x_tilde);
            Zr = Xr_tilde.*conj(Cr);
            zk = ifft(Zr);
            
            Sk(:,fDk_idx) = abs(zk).^2;
        end
    else
        for fDk_idx = 1:length(fDVtr)
            for n = 1:nStep:Nk            
                % phase estimate over the accumulation interval
                theta_jk_hat = 0;
                theta_hat = 2*pi*fDVtr(fDk_idx)*tau + theta_jk_hat;
                
                % Align C with the code in the incoming data.
                Cshift = circshift(C,n-1);
    
                % Calculate Sk and store in a matrix that represents the 2D-grid
                sk = sum(x.*exp(-1i*(theta_hat)).*Cshift);
                Sk(n, fDk_idx) = abs(sk)^2;
            end
        end
    end
end