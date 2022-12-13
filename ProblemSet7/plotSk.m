function [] = plotSk(TXID, tau, result, cfgAcquisition)
    Sk = result{TXID}.Sk;
    figure();
    h = surf(cfgAcquisition.fDVtr, tau(1:size(Sk,1))*1e6, Sk(:,:));
    set(h,'LineStyle','none')
    title('2D grid - S_k')
    xlabel('$\hat{f}_{Dk} [Hz]$','Interpreter','latex')
    ylabel('$\tau_{jk} [us]$','Interpreter','latex')
end

