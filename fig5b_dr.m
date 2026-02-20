clc; close all; clear;

% Parameters
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
L=1000;

%  depth in meters
distances = [40, 60, 80, 100];

% Wind speed (fixed)
wind_speed = 10; % Wind speed in m/s

% SNR range in dB
snrdB = 0:30;
snr = 10.^(snrdB ./ 10);

% Background radiation power reaching detector
p_s = (R_d * A * delta_t_prime * lambda * delta_lambda * omega) / (4 * h * c * delta_t);

% Initialize QBER matrix
QBER_values = zeros(length(distances), length(snrdB));

% Loop over different distances
for d = 1:length(distances)
    r = distances(d);
    
    % Turbulence effect
    sigmar2 = 37.3 * K3 * ((2 * pi / lambda)^(7/6)) * (L^(11/6));
    sigma2 = exp(((0.49 * sigmar2) / ((1 + 1.11 * sigmar2^(12/5))^(7/6))) + ((0.51 * sigmar2) / ((1 + 0.69 * sigmar2^(12/5))^(5/6)))) - 1;
    
    % Generate random values of beta
    rand_vals = 2 * pi * rand(1, 1);
    beta_vals = atan(sqrt(-sigma2 * log(rand_vals))) * (pi/180); % Calculate beta
    
    % Determine theta_1 range
    theta_min = abs(beta_vals - phi);
    theta_max = beta_vals + phi;
    theta_1 = (theta_min + (theta_max - theta_min) .* rand_vals); % Random distribution
    
    % Fresnel transmission coefficients
    theta_2 = asin((1.0 * sin(theta_1) / 1.33));
    t_p = (2 * sin(theta_2) .* cos(theta_1)) ./ (sin(theta_1 + theta_2) .* cos(theta_1 - theta_2));
    t_s = (2 * sin(theta_2) .* cos(theta_1)) ./ (sin(theta_1 + theta_2));
    
    % Photon deviation angle
    theta_dev = abs((atan(t_p ./ t_s)) - (pi / 4));
    
    % Quantum efficiency based on photon deviation angle
    eta_p = (sin(theta_dev)).^2;
    
    % Modify attenuation coefficient based on wind speed
    attenuation_factor = 1.2 + (wind_speed / 1000);
    modified_x_c = x_c * attenuation_factor;
    N_modified = (u * eta / (4 * delta_t)) * exp(-modified_x_c .* r);
    
    % Compute QBER for each SNR value
    for j = 1:length(snrdB)
        snr_linear = snr(j);
        received_signal = snr_linear * (N_modified);
        
        QBER_values(d, j) = ((received_signal * eta_p) + I_dc + p_s) ./ ...
                            ((2 * received_signal) + (2 * I_dc) + (2 * p_s));
    end
end

% Define QBER thresholds
threshold_simple_attack = 0.25; % QBER <= 25%
threshold_quantum_attack = 0.10; % QBER <= 10%

% Plot QBER vs SNR
figure;
hold on;
colors = {'r', 'g', 'b', 'm'}; % Colors for different distances

% Plot QBER for each distance
for d = 1:length(distances)
    semilogy(snrdB, QBER_values(d, :), strcat(colors{d}, 'o-'), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('L = %d m', distances(d)));
end

% Add attack thresholds
plot(snrdB, threshold_simple_attack * ones(size(snrdB)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Simple Attack Threshold (25%)');
plot(snrdB, threshold_quantum_attack * ones(size(snrdB)), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Quantum Attack Threshold (10%)');

% Axis labels and legend
xlabel('Avg. SNR (dB)');
ylabel('QBER');
title('QBER vs SNR at Different Distances');
legend('show', 'Location', 'northwest');
grid on;
ylim([1e-4 0.3]); % Adjust Y-axis for better scaling
hold off;




% Compute Bandwidth (Hz)
B = (c * delta_lambda) / (lambda^2); % 156.25 GHz

% Compute Data Rate for different SNR values
data_rate = B * log2(1 + snr); % Shannon capacity in bps



% Plot Data Rate vs SNR
figure;
semilogy(snrdB, data_rate / 1e9, 'b-o', 'LineWidth', 1.5); % Convert to Gbps
xlabel('SNR (dB)');
ylabel('Data Rate (Gbps)');
title('Data Rate vs SNR');
grid on;
