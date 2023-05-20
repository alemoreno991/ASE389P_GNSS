function [] = analizeResults(TXID, SV)
    height = 800;
    width = 600;

    fDk_hat = SV{TXID}.fDk_hat;
    theta_hat = SV{TXID}.theta_hat;
    tsk_hat = SV{TXID}.tsk_hat;
    SkdB = SV{TXID}.SkdB;
    Sk = SV{TXID}.Sk;
    CN0_hat = SV{TXID}.CN0_hat;
    t = SV{TXID}.t;
    
    %------ CN0 
    f1 = figure(1), clf
    plot(t, CN0_hat, LineWidth=1)
    title(['$CN0$ (PRN', num2str(TXID), ')'] , Interpreter='latex', FontSize=24)
    xlabel('time [s]', FontSize=24)
    ylabel('$CN0$ [dB-Hz]', Interpreter='latex', FontSize=24)
    grid on;
    axis tight
    f1.Position = [100 100 height width];

    saveas(f1, ['fig/CN0_PRN', num2str(TXID), '.png'], 'png');
    
    %------ |Sk|^2 & angle(Sk) 
    figure(2), clf
    subplot(211)
        plot(t, SkdB, LineWidth=1)
        title('$|S_k|^2$', Interpreter='latex', FontSize=24)
        xlabel('time [s]', FontSize=24)
        ylabel('$|S_k|^2$ [dB]', Interpreter='latex', FontSize=24)
        grid on;
    subplot(212)
        scatter(t, rad2deg(wrapToPi(angle((Sk)))), LineWidth=2)
        title('$angle(S_k)$', Interpreter='latex', FontSize=24)
        xlabel('time [s]', FontSize=24)
        ylabel('$angle(S_k) [deg]$', Interpreter='latex', FontSize=24)
        grid on;
    
    %------ Sk 
    f3 = figure(3), clf
    
    plot(t, abs(real(Sk)), 'k', LineWidth=2)
    hold on
    plot(t, imag(Sk), color=[0.5, 0.5, 0.5], LineWidth=2)
    xlabel('time [s]', FontSize=24)
    legend('$|real(S_k)|$', '$imag(S_k)$', Interpreter='latex', FontSize=24)
    title(['$S_k$ (PRN', num2str(TXID), ')'], Interpreter='latex', FontSize=24)
    grid on;
    axis tight
    f3.Position = [100 100 height width];

    saveas(f3, ['fig/sk_PRN', num2str(TXID), '.png'], 'png');
    
    %------ Doppler 
    f4 = figure(4), clf
    plot(t, round(fDk_hat,2), LineWidth=2)
    title(['Doppler (PRN', num2str(TXID), ')'], FontSize=24)
    xlabel('time [s]', FontSize=24)
    ylabel('$\hat{f}_D$ [Hz]', Interpreter='latex', FontSize=24)
    grid on;
    axis tight
    f4.Position = [100 100 height width];

    saveas(f4, ['fig/doppler_PRN', num2str(TXID), '.png'], 'png');
    
    %------ Phase delay
    figure(5), clf
    plot(t, unwrap(theta_hat), LineWidth=2)
    title('Phase delay', FontSize=24)
    xlabel('time [s]', FontSize=24)
    ylabel('$\hat{\theta}$ [deg]', Interpreter='latex', FontSize=24)
    grid on;
    
    %------ Code delay 
    figure(6), clf
    plot(t, round(tsk_hat*1e6,2), LineWidth=2)
    title('Code delay')
    xlabel('time [s]', FontSize=24)
    ylabel('$\hat{t}_{sk}$ [us]', Interpreter='latex', FontSize=24)
    grid on;
end

