function [] = analizeResults(TXID, SV)
    fDk_hat = SV{TXID}.fDk_hat;
    theta_hat = SV{TXID}.theta_hat;
    tsk_hat = SV{TXID}.tsk_hat;
    SkdB = SV{TXID}.SkdB;
    Sk = SV{TXID}.Sk;
    CN0_hat = SV{TXID}.CN0_hat;
    t = SV{TXID}.t;
    
    %------ CN0 
    figure(1), clf
    plot(t, CN0_hat, LineWidth=1)
    title('$CN0$', Interpreter='latex')
    xlabel('time [s]')
    ylabel('$CN0$ [dB-Hz]', Interpreter='latex')
    grid on;
    
    %------ |Sk|^2 & angle(Sk) 
    figure(2), clf
    subplot(211)
        plot(t, SkdB, LineWidth=1)
        title('$|S_k|^2$', Interpreter='latex')
        xlabel('time [s]')
        ylabel('$|S_k|^2$ [dB]', Interpreter='latex')
        grid on;
    subplot(212)
        scatter(t, rad2deg(wrapToPi(angle((Sk)))), LineWidth=2)
        title('$angle(S_k)$', Interpreter='latex')
        xlabel('time [s]')
        ylabel('$angle(S_k) [deg]$', Interpreter='latex')
        grid on;
    
    %------ Sk 
    figure(3), clf
    
    plot(t, abs(real(Sk)), 'k', LineWidth=2)
    hold on
    plot(t, imag(Sk), color=[0.5, 0.5, 0.5], LineWidth=2)
    xlabel('time [s]')
    legend('$|real(S_k)|$', '$imag(S_k)$', Interpreter='latex')
    grid on;
    
    
    %------ Doppler 
    figure(4), clf
    plot(t, round(fDk_hat,2), LineWidth=2)
    title('Doppler')
    xlabel('time [s]')
    ylabel('$\hat{f}_D$ [Hz]', Interpreter='latex')
    grid on;
    
    %------ Phase delay
    figure(5), clf
    plot(t, unwrap(theta_hat), LineWidth=2)
    title('Phase delay')
    xlabel('time [s]')
    ylabel('$\hat{\theta}$ [deg]', Interpreter='latex')
    grid on;
    
    %------ Code delay 
    figure(6), clf
    plot(t, round(tsk_hat*1e6,2), LineWidth=2)
    title('Code delay')
    xlabel('time [s]')
    ylabel('$\hat{t}_{sk}$ [us]', Interpreter='latex')
    grid on;
end

