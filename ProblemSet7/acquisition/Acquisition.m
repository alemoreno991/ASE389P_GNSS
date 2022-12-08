classdef Acquisition
    %ACQUISITION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Tc  % chip period
        Nc  % Number of chips within a code period
        varIQ 
    end
    
    methods
        function obj = Acquisition(tau, x, Ts, config)
            % Spreading code replica
            obj.Nc = 1023;
            obj.Tc = 1e-3/obj.Nc;
            code = generatePRN(34); % This satellite is not present => I'm analizing noise basically
        
            % Maximum number of code-delay samples to search through within 
            % a code period. 
            Ncode = floor(obj.Tc*obj.Nc/Ts);   
        
            % Setup accumulation configuration
            Ta = config.nTc*obj.Nc*obj.Tc; % Accumulation period
            Nk = floor(Ta/Ts);     % Number of samples in the accumulation interval
            code = repmat(code, [config.nTc 1] );
            
            % Perform non-coherent accumulation
            Sk = obj.nonCoherentAccum( tau, x, Ts, code, obj.Tc, Ncode, Nk, config );

            % Calculate the variance of the noiseIQ
            obj.varIQ = mean(Sk,"all") / (2*config.nAccum);
        end
        
        function [result] = acquire(obj, tau, x, Ts, TXID, config)
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
        %                   - threshold: threshold that indicates if a 
        %                                signal is detected.
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
        %                   - fDk_hat:   the estimated doppler [Hz]
        %
        %                   - tsk_hat:   the estimated code-delay [sec]
        %
        %                   - CN0:       the carrier-to-noise ratio [dB-Hz]
        %
        %==========================================================================
            % Spreading code replica
            code = generatePRN(TXID); 
        
            % Maximum number of code-delay samples to search through within 
            % a code period. 
            Ncode = floor(obj.Tc*obj.Nc/Ts);                                   
        
            % Setup accumulation configuration
            Ta = config.nTc*obj.Nc*obj.Tc; % Accumulation period
            Nk = floor(Ta/Ts);             % Number of samples in the 
                                           % accumulation interval
            code = repmat(code, [config.nTc 1] );
            
            % Perform non-coherent accumulation
            Sk = obj.nonCoherentAccum( tau, x, Ts, code, obj.Tc, Ncode, Nk, config);

            % Run the detection algorithm
            result = obj.computeDetection( tau, Sk, Ts, Ta, config );
        end

        function [result] = acquireFine(obj, tau, x, Ts, TXID, config)
            result = obj.acquire( tau, x, Ts, TXID, config );
        
            fStep = config.fDVtr(2) - config.fDVtr(1);
            fD_low = result.fDk_hat - 2*fStep;
            fD_high= result.fDk_hat + 2*fStep;
            config.fDVtr = fD_low:fStep/100:fD_high;
            result = obj.acquire( tau, x, Ts, TXID, config );
        end
    end

    methods (Access = private)

        function [Sk] = nonCoherentAccum(obj, tau, x, Ts, code, Tc, Ncode, Nk, config)
            % Calculate k accumulations
            Sk = zeros(Ncode, length(config.fDVtr), config.nAccum);
            for k = 1:config.nAccum
                % generate data for the k-th accumulation interval
                k_interval = (1:Nk) + Nk*(k-1);
                tau_jk = tau(k_interval);
                x_jk   = x(k_interval);
                
                % resample the spreading code at Ts
                C = oversampleSpreadingCode(code, Ts/Tc, 0, Nk, obj.Nc);
                
                % Calculate the Sk matrix for the k-th accumulation
                Sk(:,:,k) = obj.coherentAccum(tau_jk, x_jk, C, ...
                                 config.fDVtr, config.nStep, Ncode, config.method);
            end
            
            % Non-coherent accumulation
            Sk = sum(Sk,3);
        end

        function [Sk] = coherentAccum(obj, tau, x, C, fDVtr, nStep, Ncode, method)
            Sk = zeros(Ncode, length(fDVtr));
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
                    
                    Sk(:,fDk_idx) = abs(zk(1:Ncode)).^2;
                end
            else
                for fDk_idx = 1:length(fDVtr)
                    for j = 1:nStep:Ncode            
                        % phase estimate over the accumulation interval
                        theta_jk_hat = 0;
                        theta_hat = 2*pi*fDVtr(fDk_idx)*tau + theta_jk_hat;
                        
                        % Align C with the code in the incoming data.
                        Cshift = circshift(C,j-1);
            
                        % Calculate Sk and store in a matrix that 
                        % represents the 2D-grid
                        sk = sum(x.*exp(-1i*(theta_hat)).*Cshift);
                        Sk(j, fDk_idx) = abs(sk).^2;
                    end
                end
            end
        end

        function result = computeDetection(obj, tau, Sk, Ts, Ta, config)
            % Find the maximum Sk in the 2D-grid
            [max_Sk, max_idx] = max(Sk(:));
            [n, fDk_idx] = ind2sub(size(Sk), max_idx);

            % Estimated carrier to noise ratio 
            result.CN0 = 10*log10((max_Sk - 2*obj.varIQ)/(2*obj.varIQ*Ta));
            
            % Normalize Sk
            max_Sk = max_Sk/obj.varIQ;
            
            % Determine test statistic.
            s.C_N0dBHz     = result.CN0;       % Carrier to noise ratio
            s.N            = config.nAccum;    % Number of noncoherent accumulations.
            s.PfaAcq       = 1e-5;             % Desired chance of false acquisition.
            s.Ta           = Ta;               % Time of accumulation.
            s.fMax         = max(abs(config.fDVtr));% Maximum value of frequency.
            s.nCodeOffsets = floor(obj.Nc*obj.Tc/(Ts*config.nStep));  % Number of code offsets to test.
            s.ZMax         = max_Sk;           % Max value of test stat (sum of Sk^2).
            s.delZ         = 0.1;              % Resolution of output PDFs.
            
            [~,~,lambda0,~,~] = performAcqHypothesisCalcs(s);

            % Estimates
            result.isAcquired = max_Sk > 1.5*lambda0;
            result.fDk_hat = config.fDVtr(fDk_idx); % Hz
            result.tsk_hat = tau(n);                % sec

            result.Sk = Sk;
            result.sigmaIQ = sqrt(obj.varIQ);
        end
    end
end
