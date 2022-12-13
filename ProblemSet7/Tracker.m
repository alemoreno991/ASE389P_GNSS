classdef Tracker < handle
    %TRACKER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        bufferSk
        pll
        dll
        correlator

        sMix
        sigmaIQ
        Ts
        fIF
        Tc
        Ta
        Fs
        fc

        fD_hat
        theta_hat
        tsk_hat
    end
    
    methods
        function obj = Tracker(estimationSV, TXID, cfg)
            obj.sigmaIQ = cfg.sigmaIQ;
            obj.sMix = cfg.sMix;
            obj.Tc   = cfg.Tc;
            obj.fIF  = cfg.fIF;
            obj.Ts   = cfg.Ts;
            obj.Fs   = 1/obj.Ts;
            obj.fc   = cfg.fc;
            obj.Ta   = cfg.Ta;

            % |Sk|^2 running average buffer.
            obj.bufferSk = max(estimationSV.Sk(:))*ones(1, cfg.bufferSk.lenght); % used for moving average
            
            % Instantiate the loop filter
            [obj.pll.Ad, obj.pll.Bd, obj.pll.Cd, obj.pll.Dd, obj.pll.Bn_act] = ...
                configureLoopFilter( cfg.pll.Bn, cfg.Ta, cfg.pll.order );
            
            % Initial phase loop filter state.
            [eigenVtr, eigenVal] = eig(obj.pll.Ad);
            % If we take either eigenvector with eigenval=1 as `xk` we are 
            % guaranteed to obtain `xkp1 = xk` when the error `ek` is 0.  
            idx = find(diag(eigenVal) == 1, 1); % Let's take the first one we find (why not?)
            % Noticing the structure of the problem we can do the following (ad-hoc solution) 
            obj.pll.vk = -2*pi*estimationSV.fDk_hat; % estimated phase rate [rad/s]
            obj.pll.xk = obj.pll.vk/(obj.pll.Cd*eigenVtr(:,idx)) * eigenVtr(:,idx); 
            obj.pll.xkp1 = obj.pll.xk;
            
            % Delay loop parameters.
            obj.dll.Bn_target = cfg.dll.Bn;
            obj.dll.IsqQsqAvg = mean(obj.bufferSk);
            obj.dll.sigmaIQ   = cfg.sigmaIQ;
%             obj.dll.vp        = obj.sMix*obj.pll.vk/cfg.fc;
            obj.dll.vp        = 0;
            obj.dll.Tc        = obj.Tc;
            obj.dll.Ip        = 0; obj.dll.Qp        = 0;
            obj.dll.Ie        = 0; obj.dll.Qe        = 0;
            obj.dll.Il        = 0; obj.dll.Ql        = 0;

            % Correlator parameters and data.
            obj.correlator.txId        = TXID;      % PRN for target satellite
            obj.correlator.fIF         = obj.fIF;   % Intermediate frequency in Hz
            obj.correlator.Fs          = obj.Fs;    % sampling frequency, Hz
            obj.correlator.Ts          = obj.Ts;    % Sampling time interval in seconds
            obj.correlator.Ta          = obj.Ta;    % accumulation period
            obj.correlator.Tc          = obj.Tc;    % Nominal chip interval in seconds
            obj.correlator.Nc          = cfg.Nc;
            obj.correlator.nTc         = cfg.nTc;
            obj.correlator.eml         = obj.Tc/2;  % delay between early and late taps on code, sec
            obj.correlator.Nk          = floor(obj.Ta/obj.Ts); % Nk is the number of samples is one 1-ms accumulation.  It's ok for this number to be approximate

            % Use the acquisition estimation to initialize tracking
            obj.theta_hat = zeros(ceil(obj.correlator.Nk)+2,1); % Initialize the beat carrier phase estimate
            obj.fD_hat    = estimationSV.fDk_hat; % Initialize the doppler
            obj.tsk_hat   = 0; % mod(estimationSV.tsk_hat, 1e-3); % Initialize the code-delay
        end
        
        function [result] = update(obj, tVeck, xVeck)
            obj.correlator.Nk = length(tVeck);
            obj.theta_hat = obj.theta_hat(1:length(tVeck));

            % Perform correlations
            [early, prompt, late] = correlation(tVeck, xVeck, ...
                    obj.tsk_hat, obj.theta_hat, obj.correlator);

            % Update the moving window average of |Sk|^2.
            obj.bufferSk = [ abs(prompt.Sk).^2, obj.bufferSk(1:end-1) ];
            
            % single update step of a phase tracking loop
            obj.pll.Ip = real(prompt.Sk); obj.pll.Qp = imag(prompt.Sk);
            obj.pll.xk = obj.pll.xkp1;
            [obj.pll.xkp1, obj.pll.vk] = updatePll(obj.pll);

            % Single update step of a delay tracking loop 
            obj.dll.vp        =  obj.sMix * obj.pll.vk/(2*pi*obj.fc); % applying Sm/wc gain.
            obj.dll.IsqQsqAvg = mean(obj.bufferSk);
            obj.dll.Ie = real(early.Sk);  obj.dll.Qe = imag(early.Sk);
            obj.dll.Ip = real(prompt.Sk); obj.dll.Qp = imag(prompt.Sk);
            obj.dll.Il = real(late.Sk);   obj.dll.Ql = imag(late.Sk);
            vTotal = updateDll(obj.dll);
            
            % Results of the update
            obj.tsk_hat = obj.tsk_hat - vTotal * (obj.Ta);
            obj.tsk_hat = mod(obj.tsk_hat,1e-3);
            
            timeToInteg   = (0:obj.correlator.Nk+2)'*obj.Ts;
            obj.fD_hat    = obj.pll.vk/(2*pi);
            obj.theta_hat = mod(obj.theta_hat(end) + 2*pi*obj.fD_hat*(timeToInteg+obj.Ts), 2*pi);

            result.fD_hat       = obj.fD_hat;
            result.theta_hat    = obj.theta_hat(1);
            result.tsk_hat      = obj.tsk_hat;
            result.SkdB         = prompt.SkdB;
            result.Sk           = prompt.Sk;
            result.vTotal       = vTotal;
            result.CN0          = 10*log10((abs(prompt.Sk).^2 - 2*obj.sigmaIQ^2)/(2*obj.sigmaIQ^2*obj.Ta));
        end
    end
end