clc; close all; clear all; format long G

Fs = (40/7)*1e6; 
c = physconst("LightSpeed");

elapsed_rcv_samples = [ 
    200048937.4106447696685791
    200047382.089815676212310791
    200049447.80586522817611694
    200046426.1354338526725769
    200046329.94415956735610962
    200050768.35572600364685059
];

alignment_feature = 200048937; % External information (given in the instruction)

tS_seconds = [
    146046.705
    146046.703
    146046.704
    146046.697
    146046.703
    146046.710    
];

% Let compute the time it might take the signal to propagate from SV to RX
TOF = 20*1e6 / c; 
offset = tS_seconds + TOF;

% Let's use the alignment feature prime our samples and then apply the
% offset due to the time of flight (TOF)
tR_seconds = (elapsed_rcv_samples - alignment_feature)./Fs + offset;

% Calculate the pseudorange
pseudorange = c*(tR_seconds - tS_seconds)