classdef Tracker
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
            obj.Fs   = cfg.Fs;
            obj.fc   = cfg.fc;

            % |Sk|^2 running average buffer.
            obj.bufferSk = max(estimationSV.Sk)*ones(1, cfg.bufferSk.lenght); % used for moving average
            
            % Instantiate the loop filter
            [obj.pll.Ad, obj.pll.Bd, obj.pll.Cd, obj.pll.Dd, obj.pll.Bn_act] = ...
                configureLoopFilter( cfg.pll.Bn, cfg.Ta, cfg.pll.order );
            
            % Initial phase loop filter state.
            [eigenVtr, eigenVal] = eig(obj.pll.Ad);
            % If we take either eigenvector with eigenval=1 as `xk` we are 
            % guaranteed to obtain `xkp1 = xk` when the error `ek` is 0.  
            idx = find(diag(eigenVal) == 1, 1); % Let's take the first one we find (why not?)
            % Noticing the structure of the problem we can do the following (ad-hoc solution) 
            obj.pll.vk = 2*pi*estimationSV.fD; % estimated phase rate [rad/s]
            obj.pll.xk = obj.pll.vk/(s.Cd*eigenVtr(:,idx)) * eigenVtr(:,idx); 
            obj.pll.xkp1 = obj.pll.xk;
            
            % Delay loop parameters.
            obj.dll.Bn_target = cfg.dll.Bn;
            obj.dll.IsqQsqAvg = mean(obj.bufferSk);
            obj.dll.sigmaIQ   = cfg.sigmaIQ;
            obj.dll.vp        = obj.sMix*vk/cfg.fc;
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
            obj.correlator.eml         = obj.Tc/2;  % delay between early and late taps on code, sec
            obj.correlator.Nk          = floor(obj.Ta/obj.Ts); % Nk is the number of samples is one 1-ms accumulation.  It's ok for this number to be approximate

            % Use the acquisition estimation to initialize tracking
            obj.theta_hat = 0; % Initialize the beat carrier phase estimate
            obj.fD_hat    = estimationSV.fD; % Initialize the doppler
            obj.tsk_hat   = estimationSV.tsk_hat; % Initialize the code-delay
        end
        
        function [result] = update(obj, xVeck)
            % Perform correlations
            [early, prompt, late] = correlation(xVeck, obj.tsk_hat, ...
                                           obj.theta_hat, obj.fD_hat, cfg);

            % Update the moving window average of |Sk|^2.
            obj.bufferSk = [ abs(prompt.Sk).^2, obj.bufferSk(1:end-1) ];
            
            % single update step of a phase tracking loop
            obj.pll.Ip = real(prompt.Sk); obj.pll.Qp = imag(prompt.Sk);
            obj.pll.xk = obj.pll.xkp1;
            [obj.pll.xkp1, obj.pll.vk] = updatePll(obj.pll);
            obj.fD_hat   = obj.sMix * obj.pll.vk;                          % TODO: I'm not sure if the `high/low-side` mixing needs to be considered here 
            obj.theta_hat= obj.theta_hat + obj.fD_hat*obj.Ta;

            % Single update step of a delay tracking loop 
            obj.dll.vp        =  obj.sMix * obj.pll.vk/(2*pi*obj.fc); % applying Sm/wc gain.
            obj.dll.IsqQsqAvg = mean(obj.bufferSk);
            dlls.Ie = real(early.Sk);  dlls.Qe = imag(early.Sk);
            dlls.Ip = real(prompt.Sk); dlls.Qp = imag(prompt.Sk);
            dlls.Il = real(late.Sk);   dlls.Ql = imag(late.Sk);
            vTotal = updateDll(dlls);
            obj.tsk_hat = obj.tsk_hat + (1-vTotal)*(obj.Tc*obj.Nc);

            % Results of the update
            result.fD       = obj.fD;
            result.theta    = obj.theta;
            result.tsk      = obj.tsk;
        end
    end
end