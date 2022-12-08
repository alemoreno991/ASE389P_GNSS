function [early, prompt, late] = correlation(xVeck, tsk_hat, thetaHat, vk, cfg)
    % Time vector covering the accumulation
    tVec = [0:cfg.Nk-1]'*cfg.Ts;

    % Generate the phase argument of the local carrier replica
    ThetaVec = [2*pi*(cfg.fIF + vk)*tVec + thetaHat];
    
    % Generate the local carrier replica
    carrierVec = exp(-1i*ThetaVec);
    
    % Generate the +/-1-valued code (not yet oversampled)
    cacode = 2*cacodegn(cfg.txId) - 1;     
    
    % Oversample the code
    Neml = floor(cfg.eml/cfg.Fs);   % number of samples to advance/delay 
    Ndelay = floor(tsk_hat/cfg.Fs); % code-delay number of samples 
    cacode_oversampled_prompt = oversampleSpreadingCode(cacode,cfg.Ts/cfg.Tc, 0, cfg.Nk, 1023);
    cacode_oversampled_prompt = circshift(cacode_oversampled_prompt, Ndelay);
    cacode_oversampled_early  = circshift(cacode_oversampled_prompt, -Neml);
    cacode_oversampled_late   = circshift(cacode_oversampled_prompt,  Neml);
    
    % Generate the full local replica, with both code and carrier
    lpVeck = carrierVec.*cacode_oversampled_prompt;
    leVeck = carrierVec.*cacode_oversampled_early;
    llVeck = carrierVec.*cacode_oversampled_late;
    
    % Perform correlation and accumulation
    prompt.Sk = sum(xVeck.*lpVeck);
    early.Sk  = sum(xVeck.*leVeck);
    late.Sk   = sum(xVeck.*llVeck);

    % Examine the squared magnitude of Sk in dB.
    prompt.SkdB = 10*log10(abs(prompt.Sk)^2);
    early.SkdB  = 10*log10(abs(early.Sk )^2);
    late.SkdB   = 10*log10(abs(late.Sk  )^2);
end
