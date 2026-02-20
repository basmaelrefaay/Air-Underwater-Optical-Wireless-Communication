clc; close all; clear;

% --- Parameters ---
K3 = 10^-14; % Turbulence strength
lambda = 480e-9; % Wavelength in meters

R_d = 0.0125;
A = 30 * 10^-4; % Convert from cm^2 to m^2
delta_t_prime = 200 * 10^-15; % Convert ps to s
delta_t = 35 * 10^-9; % Convert ns to s
delta_lambda = 0.12 * 10^-9; % Convert nm to m
phi = 26 * (pi/180); % Convert degrees to radians
omega = 2*pi*(1 - cos(phi));
h = 6.626 * (10^-34); % Planck's constant
c = 3 * 10^8; % Speed of light
I_dc = 30; % Dark current in Hz    
u = 0.5; % Photon no. per pulse
eta = 0.3; % Quantum efficiency
x_c = 0.075; % Attenuation coefficient (m^-1)
L = 1000;

% Distances (meters)
distances = [40, 60, 80, 100];

% Wind speed (fixed)
wind_speed = 10; % m/s

% SNR range in dB
snrdB = 0:30;
snr = 10.^(snrdB ./ 10);

% Background radiation power reaching detector
p_s = (R_d * A * delta_t_prime * lambda * delta_lambda * omega) / (4 * h * c * delta_t);

% Preallocate QBER matrices
QBER_simulation = zeros(length(distances), length(snrdB));
QBER_analytic = zeros(length(distances), length(snrdB));

% --- Simulation with randomness (original) ---
for d = 1:length(distances)
    r = distances(d);
    
    sigmar2 = 37.3 * K3 * ((2 * pi / lambda)^(7/6)) * (L^(11/6));
    sigma2 = exp(((0.49 * sigmar2) / ((1 + 1.11 * sigmar2^(12/5))^(7/6))) + ...
                 ((0.51 * sigmar2) / ((1 + 0.69 * sigmar2^(12/5))^(5/6)))) - 1;
    
    % Generate random beta and theta_1
    rand_val = 2 * pi * rand(1,1);
    beta_val = atan(sqrt(-sigma2 * log(rand_val))) * (pi/180);
    
    theta_min = abs(beta_val - phi);
    theta_max = beta_val + phi;
    theta_1 = theta_min + (theta_max - theta_min) * rand_val;
    
    theta_2 = asin((1.0 * sin(theta_1) / 1.33));
    t_p = (2 * sin(theta_2) .* cos(theta_1)) ./ (sin(theta_1 + theta_2) .* cos(theta_1 - theta_2));
    t_s = (2 * sin(theta_2) .* cos(theta_1)) ./ (sin(theta_1 + theta_2));
    
    theta_dev = abs((atan(t_p ./ t_s)) - (pi / 4));
    eta_p = (sin(theta_dev)).^2;
    
    attenuation_factor = 1.2 + (wind_speed / 1000);
    modified_x_c = x_c * attenuation_factor;
    
    N_modified = (u * eta / (4 * delta_t)) * exp(-modified_x_c * r);
    
    for j = 1:length(snrdB)
        snr_linear = snr(j);
        received_signal = snr_linear * N_modified;
        
        QBER_simulation(d,j) = ((received_signal * eta_p) + I_dc + p_s) ./ ...
                               ((2 * received_signal) + (2 * I_dc) + (2 * p_s));
    end
end

% --- Analytic model (no randomness) ---
for d = 1:length(distances)
    r = distances(d);

    sigmar2 = 37.3 * K3 * ((2 * pi / lambda)^(7/6)) * (L^(11/6));
    sigma2 = exp(((0.49 * sigmar2) / ((1 + 1.11 * sigmar2^(12/5))^(7/6))) + ...
                 ((0.51 * sigmar2) / ((1 + 0.69 * sigmar2^(12/5))^(5/6)))) - 1;
    
    % Log-normal fading parameters
    mu_ln = -0.5 * log(1 + sigma2);
    sigma_ln = sqrt(log(1 + sigma2));
    lognorm_pdf = @(g) (1 ./ (g * sigma_ln * sqrt(2*pi))) .* ...
                      exp(- (log(g) - mu_ln).^2 ./ (2 * sigma_ln^2));
    
    % Use fixed average eta_p (e.g., from literature or approximate)
    eta_p_avg = 0.02; % Chosen fixed value
    
    attenuation_factor = 1.2 + (wind_speed / 1000);
    modified_x_c = x_c * attenuation_factor;
    
    N0 = (u * eta / (4 * delta_t)) * exp(-modified_x_c * r);
    
    for j = 1:length(snrdB)
        snr_linear = snr(j);
        
        integrand = @(g) ((snr_linear * N0 .* g * eta_p_avg + I_dc + p_s) ./ ...
                          (2 * snr_linear * N0 .* g + 2 * I_dc + 2 * p_s)) .* ...
                          lognorm_pdf(g);
                      
        QBER_analytic(d,j) = integral(integrand, 0, inf, 'RelTol',1e-6,'AbsTol',1e-12);
    end
end

% --- Plot both with analytic in red, simulation in blue ---
figure;
hold on;

colors_sim = {'r', 'g', 'b', 'm'};    % Simulation: blue for all distances
colors_ana = {'r', 'g', 'b', 'm'};    % Analytic: red for all distances

for d = 1:length(distances)
    semilogy(snrdB, QBER_simulation(d,:), strcat(colors_sim{d}, 'o-'), 'LineWidth', 1.5, ...
             'DisplayName', sprintf(' L=%dm', distances(d)));
end
%for d = 1:length(distances)
   % semilogy(snrdB, QBER_analytic(d,:), strcat(colors_ana{d}, 'o'), 'LineWidth', 2, ...
           %  'DisplayName', sprintf('Analytic L=%dm', distances(d)));
%end

% Thresholds
plot(snrdB, 0.25 * ones(size(snrdB)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Simple Attack Threshold (25%)');
plot(snrdB, 0.10 * ones(size(snrdB)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Quantum Attack Threshold (10%)');

xlabel('Avg. SNR (dB)');
ylabel('QBER');
title('Simulation  vs Analytic  QBER vs SNR');
legend('Location', 'northwest');
grid on;
ylim([1e-4 0.3]);
hold off;

% --- Print comparison table ---
selected_snr_indices = [1, 11, 21, 31]; % Correspond to 0,10,20,30 dB

fprintf('\nDistance(m) | SNR(dB) | QBER_simulation | QBER_analytic\n');
fprintf('---------------------------------------------------------\n');
for d = 1:length(distances)
    for idx = selected_snr_indices
        fprintf('%10d | %7.1f | %14.5e | %12.5e\n', distances(d), snrdB(idx), QBER_simulation(d, idx), QBER_analytic(d, idx));
    end
end
