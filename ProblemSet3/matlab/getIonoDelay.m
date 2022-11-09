function [delTauG] = getIonoDelay(ionodata,fc,rRx,rSv,tGPS,model)
% getIonoDelay : Return a model-based estimate of the ionospheric delay
%                experienced by a trans-ionospheric GNSS signal as it
%                propagates from a GNSS SV to the antenna of a terrestrial
%                GNSS receiver.
%
% INPUTS
%
% ionodata ------- Structure containing a parameterization of the
%                  ionosphere that is valid at time tGPS. The structure is
%                  defined differently depending on what ionospheric model
%                  is selected:
%
%                  broadcast --- For the broadcast (Klobuchar) model, ionodata
%                                is a structure containing the following fields:
%
%                       alpha0 ... alpha3 -- power series expansion coefficients
%                                            for amplitude of ionospheric delay
%                       beta0 ... beta3 -- power series expansion coefficients
%                                          for period of ionospheric plasma density 
%                                          cycle
%
%
% Other models TBD ...
%
% fc ------------- Carrier frequency of the GNSS signal, in Hz.
%
% rRx ------------ A 3-by-1 vector representing the receiver antenna position
%                  at the time of receipt of the signal, expressed in meters
%                  in the ECEF reference frame.
%
% rSv ------------ A 3-by-1 vector representing the space vehicle antenna
%                  position at the time of transmission of the signal,
%                  expressed in meters in the ECEF reference frame.
%
% tGPS ----------- A structure containing the true GPS time of receipt of
%                  the signal. The structure has the following fields:
%                  week -- unambiguous GPS week number
%                  seconds -- seconds (including fractional seconds) of the
%                  GPS week
%
% model ---------- A string identifying the model to be used in the
%                  computation of the ionospheric delay:
%                  broadcast --- The broadcast (Klobuchar) model.
%
% Other models TBD ...
%
% OUTPUTS
%
% delTauG -------- Modeled scalar excess group ionospheric delay experienced
%                  by the transionospheric GNSS signal, in seconds.
%
%+----------------------------------------------------------------------------+
% References: For the broadcast (Klobuchar) model, see IS-GPS-200F
% pp. 128-130.
%
%+============================================================================+
wgs84 = wgs84Ellipsoid('meter');
[lat,lon,h] = ecef2geodetic(wgs84, rRx(1), rRx(2), rRx(3));
[az,elev,slantRange] = ecef2aer(rSv(1), rSv(2), rSv(3), lat, lon, h, wgs84);
lambda_u = lon/180; % user geodetic longitude (semi-circles)
phi_u    = lat/180; % user geodetic latitude (semi-circles) 
A        = az/180;  % azimuth angle between user and satellite, measured clockwise 
                % positive from the true North (semi-circles) 
E        = elev/180;% elevation angle between user and satellite (semi_circle)

% earth's  central  angle  between  the  user  position  and  the  earth  
% projection  of ionospheric intersection point (semi-circles) 
Psi = 0.0137/(E+0.11) - 0.022;

% geodetic  latitude  of  the  earth  projection  of  the  ionospheric  
% intersection  point (semi-circles) 
phi_i = phi_u + Psi * cos(A*pi);
phi_i = max(-0.416, phi_i);
phi_i = min(0.416, phi_i);

% geodetic  longitude  of  the  earth  projection  of  the  ionospheric  
% intersection  point (semi-circles) 
lambda_i = lambda_u + Psi*cos(A*pi)/cos(phi_i*pi);

% geomagnetic latitude of the earth projection of the ionospheric 
% intersection point (mean ionospheric height assumed 350 km) (semi-circles) 
phi_m = phi_i + 0.064*cos((lambda_i - 1.617)*pi);

% local time (sec) 
t = 4.32*(10^4)*lambda_i + tGPS.seconds;
while t > 86400
    t = t - 86400;
end
while t < 0
    t = t + 86400;
end

PER = ionodata.broadcast.beta0 * phi_m^0 + ...
      ionodata.broadcast.beta1 * phi_m^1 + ...
      ionodata.broadcast.beta2 * phi_m^2 + ...
      ionodata.broadcast.beta3 * phi_m^3;
PER = max(72000, PER);


AMP = ionodata.broadcast.alpha0 * phi_m^0 + ...
      ionodata.broadcast.alpha1 * phi_m^1 + ...
      ionodata.broadcast.alpha2 * phi_m^2 + ...
      ionodata.broadcast.alpha3 * phi_m^3;
AMP = max(0, AMP);

% phase (radians) 
x = 2*pi*(t - 50400)/PER;

% obliquity factor (dimensionless) 
F = 1 + 16*(0.53 - E)^3;

% Estimate the ionospheric delay
if abs(x) < 1.57
    T_iono = F*(5e-9 + AMP * (1 - (x^2)/2 + (x^4)/24));
else
    T_iono = F * 5e-9;
end

if fc == 1227.44 * 1e6
    gamma = (77/60)^2;
    T_iono = gamma * T_iono;
end

delTauG = T_iono;
end


