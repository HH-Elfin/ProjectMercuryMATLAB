% plotLunaMeanAnomalyError.m
%
% Compares the mean anomaly (MA) from JPL Horizons against the value
% produced by the simulation's Keplerian propagator over time, and plots
% the error.
%
% Requires: moonData.csv (Horizons export, same format as fitLunaPrecession)
%
% Horizons CSV column header:
%   JDTDB, Calendar Date (TDB), EC, QR, IN, OM, W, Tp, N, MA, TA, A, AD, PR

clear; clc;
fclose('all');

%% ── Simulation parameters (must match celestialBodies_sol.m / updateBodyPositions.m) ──

G       = 6.674e-11;       % m^3 kg^-1 s^-2
M_earth = 5.97217e24;      % kg
M_luna  = 7.346e22;        % kg
mu      = G * (M_earth + M_luna);

a   = 3.847415e+08;                    % m
ecc = 6.983863247087546e-2;
M0  = deg2rad(2.458154691733430e2);    % rad, mean anomaly at Unix epoch

apsidalRate = 0.1643522243;            % deg/day (from fitLunaPrecession)

n_sim         = sqrt(mu / a^3);                          % rad/s, sidereal mean motion
n_anomalistic = n_sim + deg2rad(apsidalRate / 86400);    % rad/s, anomalistic mean motion

%% ── Load Horizons data ──
fid = fopen('moonData.csv', 'r');
raw = textscan(fid, '%f %*s %f %f %f %f %f %f %f %f %f %f %*[^\n]', ...
    'Delimiter', ',', 'CollectOutput', false);
fclose(fid);

jd_horizons = raw{1};
MA_horizons = raw{9};   % deg, osculating mean anomaly (advances at anomalistic rate)

% Read TA and A from Horizons (already in your textscan, just need to extract)
% Update textscan format to also grab TA (col 11) and A (col 12):
% '%f %*s %f %f %f %f %f %f %f %f %f %f %*[^\n]'
%   JD   str  EC  QR  IN  OM  W   Tp  N   MA  TA  A

TA_horizons = raw{10};  % true anomaly (deg)
A_horizons  = raw{11};  % semi-major axis (km)

% Radial distance from Horizons elements at each timestep
ecc_h = raw{2};


jd_epoch = 2440587.5;
t_sec    = (jd_horizons - jd_epoch) * 86400;
t_years  = t_sec / (365.25 * 86400);

% Kepler's equation for E at each timestep (vectorised Newton-Raphson)
M_sim = mod(M0 + n_sim .* t_sec, 2*pi);
E_sim = M_sim;
for iter = 1:100
    dE    = (M_sim - E_sim + ecc .* sin(E_sim)) ./ (1 - ecc .* cos(E_sim));
    E_sim = E_sim + dE;
    if max(abs(dE)) < 1e-10, break; end
end

r_sim      = a .* (1 - ecc .* cos(E_sim));                                       % m
r_horizons = A_horizons .* 1e3 .* (1 - ecc_h.^2) ./ (1 + ecc_h .* cosd(TA_horizons)); % m
r_err_km   = (r_sim - r_horizons) ./ 1e3;                                        % km

r_horizons = A_horizons .* 1000 .* (1 - ecc_h.^2) ./ (1 + ecc_h .* cosd(TA_horizons));  % m

% Radial distance from sim
nu_sim = 2 * atan2(sqrt(1+ecc) .* sin(E_sim/2), sqrt(1-ecc) .* cos(E_sim/2));
r_sim  = a * (1 - ecc * cos(E_sim));

%% ── Compute circular errors, then unwrap ──
% Horizons osculating MA advances at the anomalistic rate (measured from
% the instantaneously precessing periapsis), so compare against both to
% see which the data actually matches.

MA_sim_sidereal    = rad2deg(mod(M0 + n_sim         .* t_sec, 2*pi));
MA_sim_anomalistic = rad2deg(mod(M0 + n_anomalistic .* t_sec, 2*pi));

err_sidereal    = unwrapDeg(mod(MA_horizons - MA_sim_sidereal    + 180, 360) - 180);
err_anomalistic = unwrapDeg(mod(MA_horizons - MA_sim_anomalistic + 180, 360) - 180);

%% ── Report ──
p_sid  = polyfit(t_years, err_sidereal,    1);
p_anom = polyfit(t_years, err_anomalistic, 1);

fprintf('=== Luna mean anomaly error: sim vs Horizons ===\n\n');
fprintf('n_sidereal    = %.8f deg/day  (360 / 27.321582)\n', rad2deg(n_sim) * 86400);
fprintf('n_anomalistic = %.8f deg/day  (n_sid + omega_dot)\n\n', rad2deg(n_anomalistic) * 86400);
fprintf('Time range: %.2f to %.2f years\n\n', t_years(1), t_years(end));

fprintf('Sidereal comparison:\n');
fprintf('  Error at t=0:       %+.4f deg\n',   err_sidereal(1));
fprintf('  Error at t=end:     %+.4f deg\n',   err_sidereal(end));
fprintf('  Linear drift rate:  %.4f deg/year\n\n', p_sid(1));

fprintf('Anomalistic comparison:\n');
fprintf('  Error at t=0:       %+.4f deg\n',   err_anomalistic(1));
fprintf('  Error at t=end:     %+.4f deg\n',   err_anomalistic(end));
fprintf('  Linear drift rate:  %.4f deg/year\n\n', p_anom(1));

%% ── Plot ──
figure(1); clf;
tiledlayout(3, 1, 'TileSpacing', 'compact');

nexttile;
plot(t_years, err_sidereal,    'b-', 'LineWidth', 0.8); hold on;
plot(t_years, err_anomalistic, 'm-', 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 1);
xlabel('Years since Unix epoch');
ylabel('Error (deg)');
title('Sidereal vs anomalistic n — full baseline');
legend(sprintf('Sidereal  (%.3f °/yr)', p_sid(1)), ...
       sprintf('Anomalistic  (%.3f °/yr)', p_anom(1)), ...
       'Location', 'best');
grid on;

nexttile;
mask = t_years <= 2;
plot(t_years(mask), err_sidereal(mask),    'b-', 'LineWidth', 0.8); hold on;
plot(t_years(mask), err_anomalistic(mask), 'm-', 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 1);
xlabel('Years since Unix epoch');
ylabel('Error (deg)');
title('First 2 years — perturbation structure');
legend('Sidereal', 'Anomalistic', 'Location', 'best');
grid on;

nexttile;
% Whichever has lower drift — zoom on residual oscillation
[~, better] = min([abs(p_sid(1)), abs(p_anom(1))]);
if better == 1
    err_best = err_sidereal;    label = 'Sidereal';
else
    err_best = err_anomalistic; label = 'Anomalistic';
end
plot(t_years(mask), err_best(mask), 'b-', 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 1);
xlabel('Years since Unix epoch');
ylabel('Error (deg)');
title(sprintf('First 2 years — %s (lower drift)', label));
grid on;

sgtitle('Luna MA propagation error  (positive = sim lags Horizons)');

T_sid    = 27.321582 * 86400;           % sidereal period in seconds
n_true   = 2*pi / T_sid;                % rad/s
a_sidereal = (mu / n_true^2)^(1/3);    % m

T_anom   = 27.554551 * 86400;           % anomalistic period in seconds  
n_anom   = 2*pi / T_anom;
a_anomalistic = (mu / n_anom^2)^(1/3); % m

fprintf('a_sidereal    = %.6e m\n', a_sidereal);
fprintf('a_anomalistic = %.6e m\n', a_anomalistic);
fprintf('current a     = %.6e m\n', 3.847415e8);



figure(2); clf;
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
plot(t_years, r_err_km, 'b-', 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 1);
xlabel('Years since Unix epoch');
ylabel('r_{sim} − r_{Horizons}  (km)');
title('Radial distance error — full baseline');
grid on;

nexttile;
plot(t_years(mask), r_err_km(mask), 'b-', 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 1);
xlabel('Years since Unix epoch');
ylabel('Error (km)');
title('First 2 years');
grid on;

sgtitle('Luna radial distance error: sim vs Horizons');

fprintf('Radial error (km):\n');
fprintf('  RMS:  %.1f km\n',        rms(r_err_km));
fprintf('  Max:  %.1f km\n',        max(abs(r_err_km)));
fprintf('  Bias: %+.1f km\n\n',     mean(r_err_km));

function ang_unwrap = unwrapDeg(ang)
    ang_unwrap = ang;
    offset = 0;
    for k = 2:length(ang)
        d = ang(k) + offset - ang_unwrap(k-1);
        if d > 180
            offset = offset - 360;
        elseif d < -180
            offset = offset + 360;
        end
        ang_unwrap(k) = ang(k) + offset;
    end
end