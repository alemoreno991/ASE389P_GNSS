function [fDVec,tVec] = ...
simulateTrainDoppler(fc, vTrain, t0, x0, xObs, dObs, delt, N, vs)
% simulateTrainDoppler : Simulate the train horn Doppler shift scenario.
%
% INPUTS
%
% fc ------ train horn frequency, in Hz
%
% vTrain -- constant along-track train speed, in m/s
%
% t0 ------ time at which train passed the along-track coordinate x0, in
%           seconds
%
% x0 ------ scalar along-track coordinate of train at time t0, in meters
%
% xObs ---- scalar along-track coordinate of observer, in meters
%
% dObs ---- scalar cross-track coordinate of observer, in meters (i.e.,
%           shortest distance of observer from tracks)
%
% delt ---- measurement interval, in seconds
%
% N ------- number of measurements
%
% vs ------ speed of sound, in m/s
%
%
% OUTPUTS
%
% fDVec --- N-by-1 vector of apparent Doppler frequency shift measurements as
%           sensed by observer at the time points in tVec
%
% tVec ---- N-by-1 vector of time points starting at t0 and spaced by delt
%           corresponding to the measurements in fDVec
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:
%+==============================================================================+
    tVec = zeros(N-1,1);
    fDVec = zeros(N-1,1);
    
    t = 0;
    rRX = positionRX(t, xObs, dObs);
    deltaR_prev = norm(rRX - positionTX(t, x0, vTrain));
    for h = 2:N
        t = t + delt;
        TOF = norm(rRX - positionTX(t, x0, vTrain))/vs;
        deltaR_now = norm(rRX - positionTX(t, x0, vTrain));
        v_los = (deltaR_now - deltaR_prev) / delt;
        
        beta = v_los / vs;
        fr = fc / (1 + beta);
        fDVec(h-1) = fr - fc;
        tVec(h-1) = t + TOF;
        deltaR_prev = deltaR_now;
    end
end

function [rRX] = positionRX(t, x0, y0)
    rRX = [x0, y0];
end

function [rTX] = positionTX(t, x0, v0)
    x = x0 + v0*t;
    rTX = [x, 0];
end

% function [TOF] = calculateTOF(t, rRX, x0_TX, v0_TX, vs)
%     deltaR = inf;
%     TOF = norm(rRX - positionTX(t, x0_TX, v0_TX))/vs; % estimated TOF
%     while (deltaR > 1e-3)
%         deltaR= vs*TOF - norm(rRX - positionTX(t-TOF, x0_TX, v0_TX));
%         TOF = TOF + deltaR/vs;
%     end
% end