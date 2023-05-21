function [Ad,Bd,Cd,Dd,Bn_act] = configureLoopFilter(Bn_target,Ta,loopOrder)
% configureLoopFilter : Configure a discrete-time loop filter for a feedback
%                       tracking loop.
%
%
% INPUTS
%
% Bn_target ----- Target loop noise bandwidth of the closed-loop system, in
%                 Hz.
%
% Ta ------------ Accumulation interval, in seconds. This is also the loop
%                 update (discretization) interval.
%
% loopOrder ----- The order of the closed-loop system. Possible choices
%                 are 1, 2, or 3.
%
%
% OUTPUTS
%
% Ad,Bd,Cd,Dd --- Discrete-time state-space model of the loop filter.
%
% Bn_act -------- The actual loop noise bandwidth (in Hz) of the closed-loop
%                 tracking loop as determined by taking into account the
%                 discretized loop filter, the implicit integration of the
%                 carrier phase estimate, and the length of the accumulation
%                 interval.
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+

if Bn_target*Ta > 0.1
    display("Bn_target*Ta seems too big (possible discretization error). " + ...
        "I'll compute it anyways.")
end


% Loop filter transfer function
if loopOrder == 1
    % First order
    % Tune such that the loop noise bandwidth is Bn
    K = 4*Bn_target;
    num = K;
    den = 1;
    Ds = tf(num,den);
    Dz = c2d(Ds, Ta, 'zoh');    
elseif loopOrder == 2
    % Second order
    K = (8/3)*Bn_target;
    a = K/2;
    % Generate the system
    num = K*[1, a];
    den = [1, 0];
    Ds = tf(num,den);
    Dz = c2d(Ds, Ta, 'zoh');
elseif loopOrder == 3
    % Third order
    a = 1.2*Bn_target;
    b = a^2/2;
    K = 2*a;
    % Generate the system
    num = K*[1, a, b];
    den = [1, 0, 0];
    Ds = tf(num,den);
    Dz = c2d(Ds, Ta, 'zoh');
end

% Convert the loop filter to a discrete-time state-space model
[Ad,Bd,Cd,Dd]=ssdata(Dz);

% Compute the close-loop system
NCO = tf([Ta],[1 -1],Ta);       % zoh-discretized NCO
PD = 1/2*tf([1 1],[1 0],Ta);    % Phase Detector
sysOpenLoop = PD * Dz * NCO;
Hz = feedback(sysOpenLoop, 1);

% Calculate Bn_act
walias = pi/Ta;
wvec = [0:10000]'*(walias/10000);
[magvec,~] = bode(Hz, wvec);
magvec = magvec(:);
Bn_act = sum(magvec.^2)*mean(diff(wvec))/(2*pi*(magvec(1,1)^2));