clc
close all
clear all

rho = 4;
theta = 0;

% theoretical distribution
thetaML = arctanPDF(4,4,1,1);
theta_random_var = -pi/2:0.001:pi/2;
for ii = 1:length(theta_random_var)
    mu1 = rho*sin(theta);
    mu2 = rho*cos(theta);
    pdf(ii) = thetaML.calculate(mu1, mu2, theta_random_var(ii));
end

% Aproximation
for ii = 1:1e6
    n = randn() + 1i*randn();
    S(ii) = rho*exp(1i*theta) + n;
    thetaML_hat(ii) = atan(imag(S(ii)) / real(S(ii)));
end
CRLB = 1/rho;

% plots
figure()
histogram(thetaML_hat, 'Normalization', 'pdf')
hold on
plot(theta_random_var, pdf, 'b', LineWidth=3)
plot(CRLB*ones(100,1),linspace(0,1.8,100), 'r', LineWidth=3)
plot(std(thetaML_hat)*ones(100,1),linspace(0,1.8,100), 'g', LineWidth=3)
legend('$\hat{\theta}_{ML}$', '$p(\theta_{ML})$', ...
        'Cramer-Rao lower bound', '$var(\hat{\theta}_{ML})$', ...
        'Interpreter', 'latex', fontsize=24)
xlabel('$pi [rad]$','Interpreter','latex', fontsize=24)
