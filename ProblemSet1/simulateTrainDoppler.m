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
    tVec = zeros(N,1);
    fDVec = zeros(N,1);
    
    t = t0;
    for h = 1:N
        % Find rRX(t) and vRX(t)
        rRX = positionRX(t, t0, xObs, dObs);
        vRX = velocityRX(t, t0);
  
        % Calculate TOF
        TOF = calculateTOF(t, rRX, t0, x0, vTrain, vs);
        
        % Find rTX(t-TOF) and vTX(t-TOF)
        rTX = positionTX(t-TOF, t0, x0, vTrain);
        vTX = velocityTX(t-TOF, t0, vTrain);

        % Find the line-of-sight
        los = (rRX - rTX) / norm(rRX - rTX);

        % Project the velocity of TX to the line-of-sight
        vTX_los =  dot(vTX, los);
        % Project the velocity of RX to the line-of-sight
        vRX_los =  dot(vRX, los);
        % Calculate the line-of-sight velocity 
        v_los = vRX_los - vTX_los;

        % Calculate the apparent Doppler frequency shift as sensed by the
        % observer
        beta = v_los / vs;
        fr = fc / (1 + beta);
        fDVec(h) = fr - fc;
        tVec(h) = t;
        t = t + delt;
    end
end

function [rRX] = positionRX(t, t0, x0, y0)
    rRX = [x0, y0];
end

function [vRX] = velocityRX(t, t0)
    vRX = [0, 0];
end

function [rTX] = positionTX(t, t0, x0, v0)
    x = x0 + v0*(t-t0);
    rTX = [x, 0];
end

function [vTX] = velocityTX(t, t0, v0)
    vTX = [v0, 0];
end

function [TOF] = calculateTOF(t, rRX, t0, x0, vTrain, vs)
    % Initialize with a guess
    TOF = norm(rRX - positionTX(t, t0, x0, vTrain))/vs; 
    error = vs*TOF - norm(rRX - positionTX(t-TOF, t0, x0, vTrain));
    
    % Iterate until convergence
    while ( abs(error) > 1e-3 )
        TOF = TOF - error/vs;
        error = vs*TOF - norm(rRX - positionTX(t-TOF, t0, x0, vTrain));
    end
end