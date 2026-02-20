clc; close all; clear;

%% ========================== Inputs / Params ==========================
% --- Grid / channel (kept close to your Code 1) ---
K3        = 1e-13;            % underwater turbulence strength
M         = 1e5;              % Monte Carlo samples for fading
lambda_q  = 490e-9;           % wavelength (m)
distances = [40 60 80 100];   % underwater path lengths (m)
snr_dB    = 0:30;             % SNR axis (dB)
snr_lin   = 10.^(snr_dB/10);

% --- QKD / detector / environment ---
h = 6.626e-34;  c = 3e8;
A_ap = 30e-4;                 % m^2 (30 cm^2)
Delta_t      = 35e-9;         % bit period (s)
Delta_t_gate = 200e-12;       % detector gate (s)
delta_lambda = 0.12e-9;       % m
phiFOV = deg2rad(26);         % receiver FOV half-angle
Omega  = 2*pi*(1 - cos(phiFOV));
R_d    = 0.0125;              % irradiance coeff
I_dc   = 30;                  % dark counts (Hz)
eta    = 0.40;                % detector QE
u      = 0.10;                % photons per pulse (single-photon-ish)
x_c    = 0.075;               % underwater attenuation coef (1/m)

% Background photons per bit (constant)
p_s = (R_d * A_ap * Delta_t_gate * lambda_q * delta_lambda * Omega) / (4*h*c*Delta_t);

%% ===================== HOP 1: FSO (air, shore -> surface) ====================
L_air = 1000;                 %  height above sea surface (m)
DT = 5e-3;                    % Tx aperture (m)
DR = 5e-3;                    % Rx entrance pupil at the surface node (m)
theta_div = 0.5e-3;           % beam divergence (rad) ~ 0.5 mrad
alpha_atm_dBkm = 0.2;         % clear air attenuation (dB/km)

% Islam-style FSO geometric + atmospheric power transmission
T_air = 10^(-(alpha_atm_dBkm*(L_air/1000))/10) * (DR^2) / (DT + theta_div*L_air)^2;

%% ===================== HOP 2: Fiber (surface relay) ==========================
L_fiber_km     = 2.0;         % fiber length (km) — adjust to your layout
alpha_fiber_dB = 0.20;        % dB/km (G.657A2/B2 typical)
T_fiber = 10^(-(alpha_fiber_dB * L_fiber_km)/10);  % power transmission

%% ===================== HOP 3: UOWC (underwater) ==============================
% No air-water interface here (conversion at surface relay), so set extra
% polarization error for underwater hop to zero (can set small >0 if desired).
p_pol_uw = 0.0;

% Photon "base" at start of UW hop (before underwater loss & fading)
N0 = (u * eta) / (4 * Delta_t) * T_air * T_fiber;   % photons/bit arriving to UW Tx side

%% ===================== QBER vs SNR (same fading style as Code 1) ============
QBER = zeros(numel(distances), numel(snr_dB));

rng(42);  % reproducible fading
for d = 1:numel(distances)
    L = distances(d);

    % Underwater log-normal turbulence stats (your formulas)
    sigmar2 = 37.3 * K3 * ((2*pi/lambda_q)^(7/6)) * (L^(11/6));
    sigmar  = sqrt(sigmar2);
    sigma2I = exp( (0.49*sigmar2)/((1 + 1.11*sigmar^(12/5))^(7/6)) + ...
                   (0.51*sigmar2)/((1 + 0.69*sigmar^(12/5))^(5/6)) ) - 1;
    sigma2x = 0.065 * sigma2I;

    % Draw 3-branch log-amplitudes; convert to intensity gains
    hL = -sigma2x + sqrt(sigma2x).*randn(3, M);
    xa = exp(2 .* hL);                      % intensity gains per branch
    I_sum_sq = sum(xa.^2, 1);               % your MISO-RC combining metric

    % Underwater exponential loss
    T_uw = exp(-x_c * L);

    % Base photons at UW receiver *without* fading (per bit)
    N_base = N0 * T_uw;                     % scalar for this L

    for i = 1:numel(snr_lin)
        s = snr_lin(i);

        % Map SNR axis to "signal photon events" using SAME fading realizations
        S = s .* (N_base .* I_sum_sq);      % signal term aligned with your BER code

        % QBER at final receiver (trusted relays, only UW hop adds channel loss)
        q_inst = (S .* p_pol_uw + I_dc + p_s) ./ (2*S + 2*I_dc + 2*p_s);
        QBER(d, i) = mean(q_inst);
    end
end

%% ================================ Plot (linear) ==============================
figure; hold on; grid on;

cols = ['r','g','b','m'];
for d = 1:numel(distances)
    plot(snr_dB, QBER(d,:), [cols(d) 'o-'], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('L = %d m', distances(d)));
end

% Security thresholds (linear scale)
if exist('yline','file') == 2
    yline(0.25,'r','Simple Attack Threshold (25%)','LineWidth',1.5);
    yline(0.10,'b--','Quantum Attack Threshold (10%)','LineWidth',1.5);
else
    xl = xlim;                       % fallback for very old MATLAB
    plot(xl,[0.25 0.25],'r','LineWidth',1.5);
    plot(xl,[0.10 0.10],'b--','LineWidth',1.5);
end

xlabel('Avg. SNR (dB)');
ylabel('QBER');
title(sprintf('Three-hop QKD (FSO?Fiber?UOWC): L_{air}=1000 m, L_{fiber}=%.1f km', L_fiber_km));
legend('show','Location','northwest');

axis([0 30 0 0.30]);            % <<< same scale as your screenshot
set(gca,'YTick',0:0.05:0.30);   % optional finer ticks
