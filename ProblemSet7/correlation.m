function [early, prompt, late] = correlation(tVeck, xVeck, tsk_hat, thetaHat, cfg)
    
    localNk = length(tVeck);

    % Generate the phase argument of the local carrier replica
    ThetaVec = 2*pi*cfg.fIF*tVeck + thetaHat;
    
    % Generate the local carrier replica
    carrierVec = exp(-1i*ThetaVec);
    
    % Generate the +/-1-valued code (not yet oversampled)
    cacode = generatePRN(cfg.txId);     
    cacode = repmat(cacode, [cfg.nTc 1] );

    % Oversample the code
    code_delay = -tsk_hat/cfg.Tc;
    eml_delay  = -cfg.eml/(2*cfg.Tc);
    cacode_oversampled_prompt = oversampleSpreadingCode(cacode, cfg.Ts/cfg.Tc, code_delay          , localNk, cfg.Nc);
    cacode_oversampled_early  = oversampleSpreadingCode(cacode, cfg.Ts/cfg.Tc, code_delay-eml_delay, localNk, cfg.Nc);
    cacode_oversampled_late   = oversampleSpreadingCode(cacode, cfg.Ts/cfg.Tc, code_delay+eml_delay, localNk, cfg.Nc);
    
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

%     % DEBUG STUFF
%     if early.SkdB > prompt.SkdB || late.SkdB > prompt.SkdB
%         display('Sk (prompt) not the biggest')
%     end

%     idx = find(tVec>tsk_hat, floor(cfg.Nk/10), 'first');
%     figure(); 
%     plot(tVec(idx)-tsk_hat, cacode_oversampled_prompt(idx)); 
%     hold on; 
%     stem(linspace(0,1e-3,1023), cacode(1:1023))
%     legend('oversampled', 'original')
end
