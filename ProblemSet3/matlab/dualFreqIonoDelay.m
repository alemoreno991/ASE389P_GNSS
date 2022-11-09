function [t, I_L1_code, I_L2_code, I_L1_carrier, I_L2_carrier] = ...
                                                dualFreqIonoDelay(M, TXID)
week2sec = 7 * 24 * 60 * 60;
c = physconst('LightSpeed');
f_L1 = 1575.42 * 1e6;
f_L2 = 1227.60 * 1e6;
lambda_L1 = c / f_L1;
lambda_L2 = c / f_L2;
N_L1 = 0;
N_L2 = 0;

GPS_L1_CA = 0;
GPS_L2_CL = 2;

% load L1
iidum = find(M(:,13) == GPS_L1_CA & M(:,14) == TXID &...
             M(:,3) ~= 9999 & M(:,10) == 1 & M(:,11) == 0);
ORT_week = M(iidum,3);
ORT_sec_week = M(iidum,4);
ORT_frac_sec = M(iidum,5);
t_L1 = week2sec*ORT_week + ORT_sec_week + ORT_frac_sec;
rho_L1 = M(iidum,8);
phi_L1 = M(iidum,7);

% load L2
iidum = find(M(:,13) == GPS_L2_CL & M(:,14) == TXID  &...
             M(:,3) ~= 9999 & M(:,10) == 1 & M(:,11) == 0);
ORT_week = M(iidum,3);
ORT_sec_week = M(iidum,4);
ORT_frac_sec = M(iidum,5);
t_L2 = week2sec*ORT_week + ORT_sec_week + ORT_frac_sec;
rho_L2 = M(iidum,8);
phi_L2 = M(iidum,7);

% Preprocess (IDK if the L1 and L2 samples coincide temporarly)
t0 = max(t_L1(1), t_L2(1));
tf = min(t_L1(end), t_L2(end));
t = linspace(t0, tf, 1e5);
rho_L1_intrp = spline(t_L1, rho_L1, t);
phi_L1_intrp = spline(t_L1, phi_L1, t);
rho_L2_intrp = spline(t_L2, rho_L2, t);
phi_L2_intrp = spline(t_L2, phi_L2, t);

% Ionospheric delay [Code]
I_L1_code = f_L2^2 / (f_L1^2 - f_L2^2) * (rho_L2_intrp - rho_L1_intrp);
I_L2_code = f_L1^2 / (f_L2^2 - f_L1^2) * (rho_L1_intrp - rho_L2_intrp);

% Ionospheric delay [Carrier]
I_L1_carrier = f_L2^2 / (f_L1^2 - f_L2^2) * ...
       (lambda_L1*(phi_L1_intrp - N_L1) - lambda_L2*(phi_L2_intrp - N_L2));
I_L2_carrier = f_L1^2 / (f_L2^2 - f_L1^2) * ...
       (lambda_L2*(phi_L2_intrp - N_L2) - lambda_L1*(phi_L1_intrp - N_L1));

end

