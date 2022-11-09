% topPerformAcqHypothesisCalcs
%
% Top-level script for performing acquisition calculations

%----- Setup
clear; clc;
s.C_N0dBHz = 30;    
s.N = 99;            
s.PfaAcq = 0.0001;   
s.Ta = 0.001;     
s.fMax = 7000;
s.nCodeOffsets = 1023*5; 
s.ZMax = 1000;
s.delZ = 0.1;

[pZ_H0,pZ_H1,lambda0,Pd,ZVec] = performAcqHypothesisCalcs(s);

%----- Visualize the results
figure(2);
[pmax,iimax] = max(pZ_H1);
Zmax = ZVec(iimax);
clf;
ash = area(ZVec,pZ_H0);
set(get(ash,'children'), 'facecolor', 'g', 'linewidth', 2, 'facealpha', 0.5);
hold on;
ash = area(ZVec,pZ_H1);
set(get(ash,'children'), 'facecolor', 'b', 'linewidth', 2, 'facealpha', 0.5);
linemax = 1/5*max([pZ_H0;pZ_H1]);
line([lambda0,lambda0],[0,linemax], 'linewidth', 2, 'color', 'r');
xlim([0 max(Zmax*2,lambda0*1.5)]);
ylabel('Probability density');
xlabel('Z');
fs = 12;
title('GNSS Acquisition Hypothesis Testing Problem');
disp(['Probability of acquisition false alarm (PfaAcq): ' ...
      num2str(s.PfaAcq)]);
disp(['Probability of detection (Pd): ' num2str(Pd)]);
text(lambda0,linemax*1.1, ['\lambda_0 = ' num2str(lambda0) ], ...
     'fontsize',fs);
shg

%%
clc; 
%----- Execute a)
C_N0 = [27, 30, 33, 36, 39, 42, 45, 48, 51];
for idx = 1:length(C_N0)
    s.C_N0dBHz = C_N0(idx);
    flag = true;
    for N = 1:400
        s.N = N;            
        [pZ_H0,pZ_H1,lambda0,Pd,ZVec] = performAcqHypothesisCalcs(s);
        PdVec1(idx, N) = Pd;
        if (Pd > 0.95 && flag)
            N_req(idx) = N;
            flag = false;
        end
    end
end

figure()
subplot(2,1,1)
for idx = 5:length(C_N0)
    plot(PdVec1(idx,:), 'DisplayName', num2str(C_N0(idx)) )
    hold on
end
plot( 0.95*ones(size(PdVec1, 2),1), 'LineWidth', 2, 'DisplayName', 'threshold')
legend()
title('Stronger signals')
xlabel('N')
ylabel('Pd')

subplot(2,1,2)
for idx = 1:4
    plot(PdVec1(idx,:), 'DisplayName', num2str(C_N0(idx)) )
    hold on
end
plot( 0.95*ones(size(PdVec1, 2),1), 'LineWidth', 2, 'DisplayName', 'threshold')
legend()
title('Weaker signals')
xlabel('N')
ylabel('Pd')

% N_req(find(C_N0==30))

%----- Execute b)
s.N = 1;
C_N0 = 5:5:50;
TaVec = 0.001:0.001:10;
for idx = 1:length(C_N0)
    s.C_N0dBHz = C_N0(idx);
    flag = true;
    for iidx = 1:length(TaVec)
        s.Ta = TaVec(iidx);            
        [pZ_H0,pZ_H1,lambda0,Pd,ZVec] = performAcqHypothesisCalcs(s);
        PdVec2(idx, iidx) = Pd;
        if (Pd > 0.95 && flag)
            Ta_req(idx) = TaVec(iidx);
            flag = false;
        end
    end
end

figure()
subplot(2,1,1)
for idx = 5:length(C_N0)
    plot(TaVec*1e3, PdVec2(idx,:), 'DisplayName', num2str(C_N0(idx)) )
    hold on
end
plot( 0.95*ones(size(PdVec2, 2),1), 'LineWidth', 2, 'DisplayName', 'threshold')
legend()
title('Stronger signals')
xlabel('Ta [ms]')
ylabel('Pd')

subplot(2,1,2)
for idx = 1:4
    plot(TaVec*1e3, PdVec2(idx,:), 'DisplayName', num2str(C_N0(idx)) )
    hold on
end
plot( 0.95*ones(size(PdVec2, 2),1), 'LineWidth', 2, 'DisplayName', 'threshold')
legend()
title('Weaker signals')
xlabel('Ta [ms]')
ylabel('Pd')