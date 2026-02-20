clc; close all; clear;

%% ---------- Shared parameters ----------
K3 = 1e-14;                 % Turbulence strength
lambda = 480e-9;            % Wavelength (m)

R_d = 0.0125;
A = 30e-4;                  % 30 cm^2 -> m^2
delta_t_prime = 200e-15;    % s
delta_t = 35e-9;            % s
delta_lambda = 0.12e-9;     % m
phi = 26*(pi/180);          % rad
omega = 2*pi*(1 - cos(phi));
h = 6.626e-34;              % J*s
c = 3e8;                    % m/s
I_dc = 30;                  % Hz
u = 0.5;                    % photons/pulse
eta = 0.3;                  % quantum efficiency
x_c = 0.075;                % m^-1 (base attenuation)

alpha = 40*(pi/180);        % rad
n1 = 1.0;                   % air
n2 = 1.33;                  % water
L  = 1000;                  % air path (m)

% Transmission depth & wind speeds
r = [1 2 4 5 10 20 30 40 50 60 80 100 150]; % m
v = [5 10 20 40];                         % m/s

%% ---------- Background radiation ----------
p_s = (R_d * A * delta_t_prime * lambda * delta_lambda * omega) ...
      / (4 * h * c * delta_t);

%% ---------- Photon factor ----------
A_photon = (u * eta) / (4 * delta_t);

%% ---------- Turbulence statistics ----------
sigmar2 = 37.3 * K3 * ((2*pi/lambda)^(7/6)) * (L^(11/6));
sigma2  = exp(((0.49*sigmar2)/((1 + 1.11*sigmar2^(12/5))^(7/6))) + ...
              ((0.51*sigmar2)/((1 + 0.69*sigmar2^(12/5))^(5/6)))) - 1;

%% ---------- Allocate ----------
QBER_sim = zeros(length(v), length(r));
QBER_ana = zeros(length(v), length(r));

%% ---------- Simulation ----------
for i = 1:length(v)

    wind_speed = v(i);

    % Wind-modified attenuation
    x_c_mod = x_c * (1.2 + wind_speed/1000);

    % Photon counts
    N_mod  = A_photon .* exp(-x_c_mod .* r);
    N_base = A_photon .* exp(-x_c .* r);

    % Turbulence-induced angle
    u_rand = rand(1,length(r));
    beta   = atan(sqrt(-sigma2 .* log(u_rand)));

    theta_min = abs(beta - alpha);
    theta_max = beta + alpha;
    theta_1   = theta_min + (theta_max - theta_min).*rand(1,length(r));

    % Fresnel coefficients
    theta_2 = asin((n1 .* sin(theta_1)) / n2);
    t_p = (2 .* sin(theta_2) .* cos(theta_1)) ./ ...
          (sin(theta_1 + theta_2) .* cos(theta_1 - theta_2));
    t_s = (2 .* sin(theta_2) .* cos(theta_1)) ./ ...
          (sin(theta_1 + theta_2));

    % Polarization mismatch
    theta_dev = abs(atan(t_p ./ t_s) - pi/4);
    eta_p = (sin(theta_dev)).^2;

    % QBER (simulation)
    QBER_sim(i,:) = (N_base .* eta_p + I_dc + p_s) ./ ...
                    (2*N_mod + 2*I_dc + 2*p_s);
end

%% ---------- Analytic ----------
eta_p_avg = 0.01;   % average polarization error
for i = 1:length(v)
    x_c_eff = x_c * (1.2 + v(i)/1000);
    N = A_photon .* exp(-x_c_eff .* r);

    QBER_ana(i,:) = (N .* eta_p_avg + I_dc + p_s) ./ ...
                    (2*N + 2*I_dc + 2*p_s);
end

%% ---------- Plot ----------
figure;
hold on; grid on;

% Styles
lineStyles = {'k-','k-.','k:','k--'};     % Simulation
markers    = {'ko','ks','kd','k^'};       % Analytic

% Simulation (lines)
for i = 1:length(v)
    plot(r, QBER_sim(i,:), lineStyles{i}, ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('v = %d m/s (Sim)',v(i)));
end

% Analytic (markers)
for i = 1:length(v)
    plot(r, QBER_ana(i,:), markers{i}, ...
        'MarkerSize',6, ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('v = %d m/s (Ana)',v(i)));
end

% Thresholds
plot(r,0.25*ones(size(r)),'k','LineWidth',1.5,...
    'DisplayName','Simple Attack Threshold (25%)');
plot(r,0.10*ones(size(r)),'k','LineWidth',1.5,...
    'DisplayName','Quantum Attack Threshold (10%)');

xlabel('Transmission Depth (m)');
ylabel('QBER');
title('QBER vs Transmission Depth under Different Wind Speeds');
legend('Location','northwest');
ylim([0 0.3]);

hold off;
