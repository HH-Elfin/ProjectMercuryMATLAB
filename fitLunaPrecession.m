% fitLunaPrecession.m
%
% Fits linear precession rates for Luna's RAAN (OM), argument of
% periapsis (W), and mean anomaly (MA) from a JPL Horizons CSV export.
%
% Expected CSV column order (Horizons header):
%   JDTDB, Calendar Date (TDB), EC, QR, IN, OM, W, Tp, N, MA, TA, A, AD, PR
%   Col:  1                2    3   4   5   6   7   8  9  10  11 12  13  14

clear; clc;
fclose('all'); clear raw;

%% Load data
% Column 2 is a date string ("A.D. 1970-Jan-01 00:00:00.0000") which spans
% multiple comma-delimited tokens, so readmatrix cannot parse these rows.
% textscan with an explicit format skips it cleanly.
%
% Horizons CSV column layout (post-header):
%   1:  JD (TDB)
%   2:  "A.D. " (literal, skipped as string)
%   3:  date string e.g. "1970-Jan-01 00:00:00.0000" (skipped)
%   4:  EC    eccentricity
%   5:  QR    periapsis distance (km)
%   6:  IN    inclination (deg)
%   7:  OM    RAAN (deg)
%   8:  W     argument of periapsis (deg)
%   9:  Tp    time of periapsis passage (JD) — skipped
%   10: N     mean motion (deg/day) — skipped
%   11: MA    mean anomaly (deg)
%   ... remaining columns ignored

fid = fopen('moonData.csv', 'r');
raw = textscan(fid, '%f %*s %f %f %f %f %f %f %f %f %*[^\n]', ...
    'Delimiter', ',', 'CollectOutput', false);
fclose(fid);

jd  = raw{1};   % Julian date
% raw{2} = EC, raw{3} = QR, raw{4} = IN
OM  = raw{5};   % RAAN (deg)
W   = raw{6};   % Arg. of periapsis (deg)
% raw{7} = Tp, raw{8} = N
MA  = raw{9};   % Mean anomaly (deg)

% Convert JD to elapsed days since Unix epoch (JD 2440587.5 = 1970-Jan-01 00:00:00 UTC)
jd_epoch = 2440587.5;
t_days = jd - jd_epoch;

%% Unwrap angles
OM_unwrap = unwrapDeg(OM);
W_unwrap  = unwrapDeg(W);
MA_unwrap = unwrapDeg(MA);

%% Diagnostics
fprintf('raw columns at row 1: ');
for k = 1:9
    fprintf('raw{%d}=%.4f  ', k, raw{k}(1));
end
fprintf('\n\n');

fprintf('--- Diagnostics ---\n');
fprintf('Number of rows loaded: %d\n', length(jd));
fprintf('JD range:     %.3f to %.3f\n', jd(1), jd(end));
fprintf('t_days range: %.3f to %.3f\n', t_days(1), t_days(end));

fprintf('\nOM raw first 5: ');    fprintf('%.4f  ', OM(1:5));         fprintf('\n');
fprintf('OM raw last 5:  ');    fprintf('%.4f  ', OM(end-4:end));    fprintf('\n');
fprintf('OM unwrap first 5: '); fprintf('%.4f  ', OM_unwrap(1:5));   fprintf('\n');
fprintf('OM unwrap last 5:  '); fprintf('%.4f  ', OM_unwrap(end-4:end)); fprintf('\n');

fprintf('\nW raw first 5: ');    fprintf('%.4f  ', W(1:5));           fprintf('\n');
fprintf('W raw last 5:  ');    fprintf('%.4f  ', W(end-4:end));      fprintf('\n');
fprintf('W unwrap first 5: '); fprintf('%.4f  ', W_unwrap(1:5));     fprintf('\n');
fprintf('W unwrap last 5:  '); fprintf('%.4f  ', W_unwrap(end-4:end)); fprintf('\n');

fprintf('\nMA raw first 5: ');    fprintf('%.4f  ', MA(1:5));          fprintf('\n');
fprintf('MA raw last 5:  ');    fprintf('%.4f  ', MA(end-4:end));     fprintf('\n');
fprintf('MA unwrap first 5: '); fprintf('%.4f  ', MA_unwrap(1:5));    fprintf('\n');
fprintf('MA unwrap last 5:  '); fprintf('%.4f  ', MA_unwrap(end-4:end)); fprintf('\n');

fprintf('\nMax single-step jump in OM (deg): %.4f\n', max(abs(diff(OM))));
fprintf('Max single-step jump in W  (deg): %.4f\n', max(abs(diff(W))));
fprintf('Max single-step jump in MA (deg): %.4f\n', max(abs(diff(MA))));

%% Linear least-squares fits
p_OM = polyfit(t_days, OM_unwrap, 1);   % [rate (deg/day), intercept]
p_W  = polyfit(t_days, W_unwrap,  1);
p_MA = polyfit(t_days, MA_unwrap, 1);

dOM_dt = p_OM(1);   % deg/day
dW_dt  = p_W(1);    % deg/day
n_fit  = p_MA(1);   % deg/day  (mean motion from fit)

OM0_fit = p_OM(2);
W0_fit  = p_W(2);
M0_fit  = p_MA(2);  % Mean anomaly at Unix epoch, consistent with n_fit

%% Residuals
OM_residual = OM_unwrap - polyval(p_OM, t_days);
W_residual  = W_unwrap  - polyval(p_W,  t_days);
MA_residual = MA_unwrap - polyval(p_MA, t_days);

OM_rms = rms(OM_residual);
W_rms  = rms(W_residual);
MA_rms = rms(MA_residual);

%% Report
n_expected = 360 / 27.321582;  % deg/day from sidereal period

fprintf('\n=== Luna element fit results ===\n\n');
fprintf('Baseline: %.1f days (%.2f years)\n\n', t_days(end) - t_days(1), (t_days(end) - t_days(1)) / 365.25);

fprintf('RAAN (OM):\n');
fprintf('  Rate (fitted)  = %+.8f deg/day\n', dOM_dt);
fprintf('  Epoch value    = %.6f deg  (Horizons raw: 344.485655 deg)\n', OM0_fit);
fprintf('  RMS residual   = %.4f deg\n\n', OM_rms);

fprintf('Arg. of periapsis (W):\n');
fprintf('  Rate (fitted)  = %+.8f deg/day\n', dW_dt);
fprintf('  Epoch value    = %.6f deg  (Horizons raw: 327.931418 deg)\n', W0_fit);
fprintf('  RMS residual   = %.4f deg\n\n', W_rms);

fprintf('Mean anomaly (MA):\n');
fprintf('  n (fitted)     = %+.8f deg/day\n', n_fit);
fprintf('  n (expected)   =  %.8f deg/day  (360 / 27.321582 days)\n', n_expected);
fprintf('  n difference   = %+.4e deg/day  (%.4f sec/orbit)\n', ...
    n_fit - n_expected, (n_fit - n_expected) / n_expected * 27.321582 * 86400);
fprintf('  M0 (fitted)    = %.6f deg  (Horizons raw: 245.815469 deg)\n', M0_fit);
fprintf('  RMS residual   = %.4f deg\n\n', MA_rms);

fprintf('--- Copy these into celestialBodies_sol.m ---\n');
fprintf('meanAnomaly           = %.10f, ... %% deg, at Unix epoch, consistent with fitted n\n', M0_fit);
fprintf('raanPrecessionRate    = %.10f, ... %% deg/day\n', dOM_dt);
fprintf('apsidalPrecessionRate = %.10f, ... %% deg/day\n', dW_dt);
fprintf('\n(Also update RAAN and W epoch values if using fit intercepts)\n');
fprintf('rightAscensionOfAscendingNode = %.6f, ... %% deg, fit intercept at Unix epoch\n', OM0_fit);
fprintf('argumentOfPeriapsis           = %.6f, ... %% deg, fit intercept at Unix epoch\n', W0_fit);

%% Plots — figure 1: OM and W (existing)
figure(1); clf;
tiledlayout(2, 2);

nexttile;
plot(t_days / 365.25, OM_unwrap, 'b.', 'MarkerSize', 2); hold on;
plot(t_days / 365.25, polyval(p_OM, t_days), 'r-', 'LineWidth', 1.5);
xlabel('Years since epoch'); ylabel('OM (deg, unwrapped)');
title(sprintf('RAAN fit  (%.4f deg/day)', dOM_dt));
legend('Data', 'Linear fit'); grid on;

nexttile;
plot(t_days / 365.25, OM_residual, 'b.', 'MarkerSize', 2);
xlabel('Years since epoch'); ylabel('Residual (deg)');
title(sprintf('RAAN residuals  RMS = %.3f deg', OM_rms));
yline(0, 'r'); grid on;

nexttile;
plot(t_days / 365.25, W_unwrap, 'b.', 'MarkerSize', 2); hold on;
plot(t_days / 365.25, polyval(p_W, t_days), 'r-', 'LineWidth', 1.5);
xlabel('Years since epoch'); ylabel('W (deg, unwrapped)');
title(sprintf('Arg. periapsis fit  (%.4f deg/day)', dW_dt));
legend('Data', 'Linear fit'); grid on;

nexttile;
plot(t_days / 365.25, W_residual, 'b.', 'MarkerSize', 2);
xlabel('Years since epoch'); ylabel('Residual (deg)');
title(sprintf('W residuals  RMS = %.3f deg', W_rms));
yline(0, 'r'); grid on;

sgtitle('Luna orbital element precession — linear fit to Horizons data');

%% Plots — figure 2: MA
figure(2); clf;
tiledlayout(1, 2);

nexttile;
plot(t_days / 365.25, MA_unwrap, 'b.', 'MarkerSize', 2); hold on;
plot(t_days / 365.25, polyval(p_MA, t_days), 'r-', 'LineWidth', 1.5);
xlabel('Years since epoch'); ylabel('MA (deg, unwrapped)');
title(sprintf('Mean anomaly fit  (n = %.6f deg/day)', n_fit));
legend('Data', 'Linear fit'); grid on;

nexttile;
plot(t_days / 365.25, MA_residual, 'b.', 'MarkerSize', 2);
xlabel('Years since epoch'); ylabel('Residual (deg)');
title(sprintf('MA residuals  RMS = %.3f deg', MA_rms));
yline(0, 'r'); grid on;

sgtitle('Luna mean anomaly — linear fit to Horizons data');

%% Helper function
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