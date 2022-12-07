classdef Acquisition
    %ACQUISITION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Tc
        fcL1 
    end
    
    methods
        function obj = Acquisition(inputArg1,inputArg2)
            obj.Tc   = 1e-3/1023; % Chipping period.
            obj.fcL1 = 1575.42e6; % L1 carrier frequency
            
            doppler_max = 10e3;          % Doppler search bounds.
            doppler_del = 1/(Taccum);    % Doppler search resolution. (Consistent with acquisition hypothesis stat calculations.)


        end
        
        function outputArg = method1(obj, Ta)
            Nk = round(Taccum*Fs);  % Number of samples in accum. interval.

        end
    end
end

