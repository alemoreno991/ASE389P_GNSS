function [Ccoh] = computeCoherence(DeltaThetaVec,N)
% computeCoherence : Compute the value of the discrete-time coherence function
%                    Ccoh(N).
%
%
% INPUTS
%
% DeltaThetaVec ----- Ns-by-1 vector representing a sampled carrier phase
%                     error time history, in rad.
%
% N ----------------- The number of samples that will be used to evaluate
%                     the coherence Ccoh(N).
%
%
% OUTPUTS
%
% Ccoh -------------- The value of the discrete-time coherence function for
%                     the first N samples of DeltaThetaVec.
%
%+------------------------------------------------------------------------------+
% References:
%
%+==============================================================================+
Ccoh = abs(sum(exp(i*DeltaThetaVec(1:N))/N));