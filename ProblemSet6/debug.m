clc; close all; clear all;

t = 0:1e-6:3;
x = simulatePhaseShift(t, 0, 'custom');

plot(t,x)

function thetaVtr = simulatePhaseShift(t, theta0, model)
% simulatePhaseShift: simulates a phase shift of the carrier 
%
%
% INPUTS
%
% t      ------------ Time vector
%
% theta0 ------------ Initial phase 
%
% model   ----------- Flag that indicates the model to implement
%
%                   Options:
%                       - 'ramp'
%                       - 'quadratic'
%
% OUTPUTS
%
% thetaVtr -------- Phase shift corresponding to each time step.
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+
% Make phase shift vector.
    switch model
        case 'custom'
            t1 = t(1:floor(length(t)/3));
            tF = t1(end);
            ramp1  = 2*pi*t1;
            t2 = t(floor(length(t)/3)+1: floor(2*length(t)/3)) - tF;
            tF = tF + t2(end);
            quadratic = ramp1(end) - pi*t2.^2;
            t3 = t(floor(2*length(t)/3)+1: end) - tF;
            ramp2  = quadratic(end) + 2*pi*t3; 
            
            thetaVtr = [ramp1, quadratic, ramp2];    
    end
    thetaVtr = thetaVtr + theta0;
end