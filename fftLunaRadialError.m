% fftLunaRadialError.m
%
% Computes and plots the FFT of the Luna radial distance error (r_sim - r_Horizons)
% to identify the dominant perturbation frequencies.
%
% Run after plotLunaMeanAnomalyError.m, or standalone (re-derives the error).
%
% Requires: moonData.csv

clear; clc;
fclose('all');

% Reset figure defaults
set(groot, 'DefaultAxesColor',          'white');
set(groot, 'DefaultAxesXColor',         'black');
set(groot, 'DefaultAxesYColor',         'black');
set(groot, 'DefaultAxesZColor',         'black');
set(groot, 'DefaultAxesGridColor',      'black');
set(groot, 'DefaultFigureColor',        'white');
set(groot, 'DefaultTextColor',          'black');
%set(groot, 'DefaultAxesFontColor',      'black');

%% ── Simulation parameters ──
G       = 6.674e-11;
M_earth = 5.97217e24;
M_luna  = 7.346e22;
mu      = G * (M_earth + M_luna);

a   = 3.847415e+08;
ecc = 6.983863247087546e-2;
M0  = deg2rad(2.458154691733430e2);

n_sim = sqrt(mu / a^3);

%% ── Load Horizons data ──
fid = fopen('moonData.csv', 'r');
raw = textscan(fid, '%f %*s %f %f %f %f %f %f %f %f %f %f %*[^\n]', ...
    'Delimiter', ',', 'CollectOutput', false);
fclose(fid);

jd_horizons = raw{1};
ecc_h       = raw{2};
TA_horizons = raw{10};
A_horizons  = raw{11};

jd_epoch = 2440587.5;
t_sec    = (jd_horizons - jd_epoch) * 86400;
t_years  = t_sec / (365.25 * 86400);

%% ── Compute radial error ──
M_sim = mod(M0 + n_sim .* t_sec, 2*pi);
E_sim = M_sim;
for iter = 1:100
    dE    = (M_sim - E_sim + ecc .* sin(E_sim)) ./ (1 - ecc .* cos(E_sim));
    E_sim = E_sim + dE;
    if max(abs(dE)) < 1e-10, break; end
end

r_sim      = a .* (1 - ecc .* cos(E_sim));
r_horizons = A_horizons .* 1e3 .* (1 - ecc_h.^2) ./ (1 + ecc_h .* cosd(TA_horizons));
r_err_km   = (r_sim - r_horizons) ./ 1e3;

%% ── Check sampling uniformity ──
dt_sec = diff(t_sec);
fprintf('Sampling interval: min=%.1f s, max=%.1f s, mean=%.1f s\n', ...
    min(dt_sec), max(dt_sec), mean(dt_sec));

% FFT requires uniform sampling. Resample onto a uniform grid if needed.
dt_uniform = mean(dt_sec);
t_uniform  = (t_sec(1) : dt_uniform : t_sec(end))';
r_err_uniform = interp1(t_sec, r_err_km, t_uniform, 'linear');

N  = length(r_err_uniform);
fs = 1 / dt_uniform;              % samples per second

%% ── FFT ──
Y    = fft(r_err_uniform);
P2   = abs(Y / N);
P1   = P2(1 : floor(N/2) + 1);
P1(2:end-1) = 2 * P1(2:end-1);   % single-sided amplitude spectrum

f_hz   = fs * (0 : floor(N/2)) / N;              % Hz
f_cpd  = f_hz * 86400;                            % cycles per day
T_days = 1 ./ f_cpd;                              % period in days

% Mask DC and unreasonably short periods (<5 days) for peak finding
valid = T_days >= 5 & T_days <= 365*20;
P1_valid = P1;
P1_valid(~valid) = 0;

% Find top 8 peaks
[pks, locs] = localMaxima(P1_valid, 8);

fprintf('\n=== Top spectral peaks ===\n');
fprintf('  Rank  Period (days)   Amplitude (km)   Known source\n');
fprintf('  -------------------------------------------------------\n');

known = {27.55,  'Anomalistic month'; ...
         27.32,  'Sidereal month'; ...
         29.53,  'Synodic month'; ...
         31.81,  'Evection (~31.8 d)'; ...
         14.77,  'Variation (~14.8 d)'; ...
         206.7,  'Annual equation (~207 d)'; ...
         365.25, 'Annual (solar)'; ...
         182.6,  'Semi-annual'; ...
         13.66,  'Fortnightly'; ...
         6798,   'Nodal (18.6 yr)'};

for k = 1:length(locs)
    T_k = T_days(locs(k));
    A_k = pks(k);
    % Find closest known period
    diffs = cellfun(@(x) abs(x - T_k), known(:,1));
    [d, ki] = min(diffs);
    if d / T_k < 0.08
        label = known{ki, 2};
    else
        label = '(unknown)';
    end
    fprintf('  %4d  %13.2f   %15.1f   %s\n', k, T_k, A_k, label);
end

%% ── Plot ──
figure(1); clf;
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
plot(T_days, P1, 'b-', 'LineWidth', 0.8);
xlim([5 365]);
xlabel('Period (days)');
ylabel('Amplitude (km)');
title('Single-sided amplitude spectrum of radial error — short periods');
grid on;

% Annotate peaks
hold on;
for k = 1:length(locs)
    T_k = T_days(locs(k));
    if T_k <= 365*3
        xline(T_k, 'r--', sprintf('%.1f d', T_k), ...
            'LabelVerticalAlignment', 'bottom', 'LineWidth', 1);
    end
end

nexttile;
plot(T_days / 365.25, P1, 'b-', 'LineWidth', 0.8);
xlim([0 20]);
xlabel('Period (years)');
ylabel('Amplitude (km)');
title('Long-period content (0–20 years)');
grid on;

hold on;
for k = 1:length(locs)
    T_k = T_days(locs(k));
    if T_k / 365.25 <= 20
        xline(T_k / 365.25, 'r--', sprintf('%.2f yr', T_k/365.25), ...
            'LabelVerticalAlignment', 'bottom', 'LineWidth', 1);
    end
end

sgtitle('FFT of Luna radial distance error (r_{sim} − r_{Horizons})');

function [pks, locs] = localMaxima(x, n)
    % Returns the n largest local maxima in x, sorted by amplitude descending
    is_peak = x(2:end-1) > x(1:end-2) & x(2:end-1) > x(3:end);
    idx = find(is_peak) + 1;  % offset for the trimmed edges
    [sorted_vals, order] = sort(x(idx), 'descend');
    take = min(n, length(order));
    pks  = sorted_vals(1:take);
    locs = idx(order(1:take));
end