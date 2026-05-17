% fitLunaPrecession.m
%
% Fits linear precession rates for Luna's RAAN (OM) and argument of
% periapsis (W) from a JPL Horizons CSV export.
%
% Expected CSV column order (no header row):
%   1: JD (TDB)
%   2: Date string (ignored)
%   3: EC   - eccentricity
%   4: QR   - periapsis distance (km)
%   5: IN   - inclination (deg)
%   6: OM   - RAAN (deg)
%   7: W    - argument of periapsis (deg)
%   ... (remaining columns ignored)

clear; clc;
fclose('all'); clear raw;

%% Load data
% Column 2 is a date string ("A.D. 1970-Jan-01 00:00:00.0000") which spans
% multiple comma-delimited tokens, so readmatrix cannot parse these rows.
% textscan with an explicit format skips it cleanly.
%
% Horizons CSV column layout:
%   1:  JD (TDB)
%   2:  "A.D. " (literal)
%   3:  date string e.g. "1970-Jan-01 00:00:00.0000"
%   4:  EC    eccentricity
%   5:  QR    periapsis distance (km)
%   6:  IN    inclination (deg)
%   7:  OM    RAAN (deg)
%   8:  W     argument of periapsis (deg)
%   ... remaining columns ignored

fid = fopen('moonData.csv', 'r');
raw = textscan(fid, '%f %*s %f %f %f %f %f %*[^\n]', ...
    'Delimiter', ',', 'CollectOutput', false);
fclose(fid);

jd = raw{1};   % Julian date
% raw{2} = EC, raw{3} = QR, raw{4} = IN
OM = raw{5};
W  = raw{6};

% Convert JD to elapsed days since Unix epoch (JD 2440587.5 = 1970-Jan-01 00:00:00 UTC)
jd_epoch = 2440587.5;
t_days = jd - jd_epoch;

%% Unwrap angles to remove 0/360 discontinuities before fitting
% OM and W are both in [0, 360); unwrapping in radians then converting back
% gives a continuous signal suitable for linear regression.
OM_unwrap = unwrapDeg(OM);
W_unwrap  = unwrapDeg(W);

fprintf('raw columns at row 1: ');
for k = 1:6
    fprintf('raw{%d}=%.4f  ', k, raw{k}(1));
end
fprintf('\n');

fprintf('--- Diagnostics ---\n');
fprintf('Number of rows loaded: %d\n', length(jd));
fprintf('JD range: %.3f to %.3f\n', jd(1), jd(end));
fprintf('t_days range: %.3f to %.3f\n', t_days(1), t_days(end));

fprintf('\nOM raw first 5: '); fprintf('%.4f  ', OM(1:5)); fprintf('\n');
fprintf('OM raw last 5:  '); fprintf('%.4f  ', OM(end-4:end)); fprintf('\n');
fprintf('OM unwrap first 5: '); fprintf('%.4f  ', OM_unwrap(1:5)); fprintf('\n');
fprintf('OM unwrap last 5:  '); fprintf('%.4f  ', OM_unwrap(end-4:end)); fprintf('\n');

fprintf('\nW raw first 5: '); fprintf('%.4f  ', W(1:5)); fprintf('\n');
fprintf('W raw last 5:  '); fprintf('%.4f  ', W(end-4:end)); fprintf('\n');
fprintf('W unwrap first 5: '); fprintf('%.4f  ', W_unwrap(1:5)); fprintf('\n');
fprintf('W unwrap last 5:  '); fprintf('%.4f  ', W_unwrap(end-4:end)); fprintf('\n');

fprintf('\nMax single-step jump in OM (deg): %.4f\n', max(abs(diff(OM))));
fprintf('Max single-step jump in W  (deg): %.4f\n', max(abs(diff(W))));

%% Linear least-squares fit: angle = angle0 + rate * t
% Using polyfit (degree 1) for simplicity.

p_OM = polyfit(t_days, OM_unwrap, 1);   % [rate (deg/day), intercept]
p_W  = polyfit(t_days, W_unwrap,  1);

dOM_dt = p_OM(1);   % deg/day
dW_dt  = p_W(1);    % deg/day

OM0_fit = p_OM(2);  % Fitted value at epoch (should be close to Horizons epoch value)
W0_fit  = p_W(2);

%% Residuals
OM_residual = OM_unwrap - polyval(p_OM, t_days);
W_residual  = W_unwrap  - polyval(p_W,  t_days);

OM_rms = rms(OM_residual);
W_rms  = rms(W_residual);

%% Report
fprintf('=== Luna precession fit results ===\n\n');
fprintf('Baseline: %.1f days (%.2f years)\n\n', t_days(end) - t_days(1), (t_days(end) - t_days(1)) / 365.25);

fprintf('RAAN (OM):\n');
fprintf('  Rate       = %+.6f deg/day\n', dOM_dt);
fprintf('  Epoch value= %.6f deg  (Horizons: 344.485655 deg)\n', OM0_fit);
fprintf('  RMS residual = %.4f deg\n\n', OM_rms);

fprintf('Arg. of periapsis (W):\n');
fprintf('  Rate       = %+.6f deg/day\n', dW_dt);
fprintf('  Epoch value= %.6f deg  (Horizons: 327.931418 deg)\n', W0_fit);
fprintf('  RMS residual = %.4f deg\n\n', W_rms);

fprintf('--- Copy these into celestialBodies_sol.m ---\n');
fprintf('raanPrecessionRate    = %.10f, ... %% deg/day\n', dOM_dt);
fprintf('apsidalPrecessionRate = %.10f, ... %% deg/day\n', dW_dt);

%% Plot
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

function ang_unwrap = unwrapDeg(ang)
    ang_unwrap = ang;
    offset = 0;
    for k = 2:length(ang)
        diff = ang(k) + offset - ang_unwrap(k-1);
        if diff > 180
            offset = offset - 360;
        elseif diff < -180
            offset = offset + 360;
        end
        ang_unwrap(k) = ang(k) + offset;
    end
end