clc; close all; clear;

% Constants
R_d = 0.0125;
A = 30 * 10^-4; % Convert from cm^2 to m^2
delta_t_prime = 200 * 10^-15; % Convert ps to s
delta_t = 35 * 10^-9; % Convert ns to s
delta_lambda = 0.12 * 10^-9; % Convert nm to m
lambda = 480 * 10^-9; % Convert nm to m
phi = 26 *(pi/180) ; % Convert degrees to radians
omega = 2*pi*(1 - cos(phi));
h = 6.626 * (10^-34); % Planck's constant
c = 3 * 10^5; % Speed of light
I_dc = 30; % Dark current in Hz    
u = 0.5; % Photon no. per pulse
eta = 0.3; % Quantum efficiency
x_c = 0.075; % Attenuation coefficient (m^-1)
r = linspace(0,150,151); % Distance range
v = [5, 10, 20, 40]; % Wind speeds
alpha = 50*(pi/180); % Angle in rad
n1 = 1.0; % Refractive index of air
n2 = 1.33; % Refractive index of water

% Calculate background radiation power reaching detector
p_s = (R_d * A * delta_t_prime * lambda * delta_lambda * omega) / (4 * h * c * delta_t);

% Function to calculate the number of photons arriving at the receiver
N = ((u * eta) .* exp(-x_c .* r))/ (4 * delta_t);

% Loop over wind speeds
for i = 1:length(v)
    wind_speed = v(i);
    
    % Calculate sigma^2 based on wind speed
    sigma2 = (0.003+ 0.000512) .* wind_speed;
    
    % Generate random values of beta based on the distribution
    rand_vals = 2*pi.*rand(1, 151); % Generate random numbers in [0, 1]
    beta_vals = atan(sqrt(-sigma2 .* log(rand_vals)))*(pi/180) ; % Calculate beta
    
    % Determine theta_1 range based on beta and alpha
    theta_min = abs(beta_vals - alpha);
    theta_max = beta_vals + alpha;
    theta_1 = (theta_min + (theta_max - theta_min) .*rand_vals); % Random distribution
    
    % Fresnel transmission coefficients
    theta_2 = (asin((n1 .* sind(theta_1) / n2)))*(pi/180);
    t_p = (2 * sin(theta_2) .* cos(theta_1)) ./ (sin(theta_1 + theta_2) .* cos(theta_1 - theta_2));
    t_s = (2 * sin(theta_2) .* cos(theta_1)) ./ (sin(theta_1 +theta_2));
    
    % Calculate photon deviation angle
    theta_dev = abs((atan(t_p ./ t_s))- (pi / 4));
    
    % Quantum efficiency based on photon deviation angle
    eta_p = (sin(theta_dev)).^2;
    
    % Modify attenuation coefficient based on wind speed
    attenuation_factor = 1.2+(wind_speed/1000 );  % Simple model of attenuation
    modified_x_c = x_c .* attenuation_factor;   % Modify attenuation coefficient based on wind speed
    N_modified =( u * eta / (4 * delta_t)) .* exp(-modified_x_c .* r); % Adjust photon count
    
    % Store QBER values for the current wind speed in QBER_values matrix
    QBER_values(i, :) =((N .* eta_p )+ I_dc + p_s) ./ ((2 .* N_modified) + (2 * I_dc) + (2 * p_s));
end

% Define QBER thresholds for safe transmission depths
threshold_simple_attack = 0.25; % QBER <= 25%
threshold_quantum_attack = 0.10; % QBER <= 10%

% Plot QBER vs Transmission Depth
figure;
hold on;
colors = {'k-', 'k-.', 'k:', 'k--'}; % Colors and line styles for different wind speeds

% Plot QBER for each wind speed with the new colors and styles
for i = 1:length(v)
    plot(r, QBER_values(i, :), colors{i}, 'LineWidth', 1.5, 'DisplayName', sprintf('v = %d m/s', v(i)));
end

% Add thresholds as horizontal lines using plot
plot(r, threshold_simple_attack * ones(size(r)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Simple Attack Threshold (25%)');
plot(r, threshold_quantum_attack * ones(size(r)), 'k.', 'LineWidth', 1.5, 'DisplayName', 'Quantum Attack Threshold (10%)');

% Axis labels and legend
xlabel('Transmission Depth (m)');
ylabel('QBER');
title('QBER vs Transmission Depth under Different Wind Speeds');
legend('show', 'Location', 'northwest');
grid on;

% Set the y-axis limits from 0.0 to 0.3
ylim([0.0 0.3]);

hold off;
