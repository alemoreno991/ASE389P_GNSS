classdef InputSignal
    %InputSignal 
    
    properties
        Tfull
        Fs
        N
        X
        bandpass
        high_side_mix
    end
    
    methods
        function [obj] = InputSignal(filename, Tfull, Fs, bandpass, mixing)
            %----- Setup
            obj.Tfull = Tfull;
            obj.Fs = Fs;
            obj.N = floor(Fs*Tfull/16)*16;
            obj.bandpass.flag = bandpass.flag;

            %----- Load data
            if strcmp(filename, 'rawintegersamples_fe.bin')
                obj.bandpass.fIF  = bandpass.fIF;
                
                stream = 1;         % Data stream (between 1 and numStreams)
                numStreams = 4;     % Number of data streams 
                tSeek = 0;          % Seek time into data (sec)

                fid = fopen(filename, 'r', 'n');
                Ns = floor(Tfull*Fs);
                % numStreams bytes per sample, one for each data stream
                seekOffset = floor(tSeek*Fs)*numStreams;
                status = fseek(fid,seekOffset,-1);
                if(status == -1)
                  error('tSeek beyond file limit');
                end
                Y = fread(fid, [numStreams,Ns], 'int8')';
                if(length(Y(:,1)) < Ns)
                  error('Insufficient data');
                end
                obj.X = Y(:,stream);
            else
                fid = fopen(filename,'r','l');
                
                % Configure Baseband/Bandpass representation
                if obj.bandpass.flag == true
                    obj.bandpass.fIF  = bandpass.fIF;
                    [X, ~] = binloadSamples(fid, obj.N, 'dual');
                    obj.X = X(:,1);
                else
                    fid = fopen('niData01head_5MHz.bin','r','l');
                    X = fread(fid, [2, obj.N], 'int16')';
                    obj.X = X(:,1) + 1j*X(:,2);
                end
            end
            fclose(fid);

            % Determinte if `high` or `low` side mixing was used
            if strcmp(mixing, 'high')
                obj.high_side_mix = true;
            else
                obj.high_side_mix = false;
            end
        end

        function [tauVec, xVec, Ts, sMix] = getBasebandRepresentation(obj)
            Ts = 1/obj.Fs;

            if obj.bandpass.flag == true
                % Convert signal to its complex baseband representation
                [IVec,QVec] = if2iq(obj.X, Ts, obj.bandpass.fIF);
                xVec = IVec - 1j*QVec;        
                Ts = Ts*2;
            else
                % complex baseband representation
                xVec = obj.X; 
            end

            tauVec = Ts*(0:1:length(xVec))';
            
            if obj.high_side_mix
                sMix = -1;
            else
                sMix = 1;
            end
        end

        function [tauVec, xVec, Ts, sMix] = getBandpassRepresentation(obj, fIF)
            Ts = 1/obj.Fs;

            if obj.bandpass.flag == false
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % TODO: not very well tested (I basically trusted my logic)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Convert signal to its bandpass representation
                Ts = 1/obj.Fs;
                xVec = iq2if(real(obj.X), imag(obj.X), Ts, fIF);
                Ts = Ts/2;
            else
                % bandpass representation
                xVec = obj.X;
            end

            tauVec = Ts*(0:1:length(xVec))';

            if obj.high_side_mix
                sMix = -1;
            else
                sMix = 1;
            end
        end
    end
end