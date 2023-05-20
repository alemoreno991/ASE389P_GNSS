% loadRawSamples
%
% Raw integer samples are read from file byte by byte as
%
% [Stream1_Integer1 Stream2_Integer1 Stream3_Integer1 ... 
%  StreamN_Integer1 Stream1_Integer2 Stream2_Integer2 ... ]
%
% Streams are ordered according to the CIRCBUFF_STREAM_IDX configuration
% for each bank; typically the ordering is
% {L1, L1_ALT1, L1_ALT2, ..., L1_ALTM, L2, L2_ALT1, ...}

clear;clc;
%----- Setup
filename = 'rawintegersamples_fe.bin';

stream = 1;         % Data stream (between 1 and numStreams)
numStreams = 4;     % Number of data streams 
Tfull = 60;          % Interval of data to load (sec)
fs = 9.6e6;      % Sampling frequency (Hz)
tSeek = 0;          % Seek time into data (sec)

%----- Load data
fid = fopen(filename, 'r', 'n');
Ns = floor(Tfull*fs);
% numStreams bytes per sample, one for each data stream
seekOffset = floor(tSeek*fs)*numStreams;
status = fseek(fid,seekOffset,-1);
if(status == -1)
  error('tSeek beyond file limit');
end
Y = fread(fid, [numStreams,Ns], 'int8')';
fclose(fid);
if(length(Y(:,1)) < Ns)
  error('Insufficient data');
end
Y = Y(:,stream);
