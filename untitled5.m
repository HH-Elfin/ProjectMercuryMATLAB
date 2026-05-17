% Known Horizons values at 2026-May-17 (your diagnostic date)
M_horizons  = deg2rad(351.64);
nu_horizons = meanToTrue(M_horizons, 0.0549); % use your Kepler solver

% What r does Horizons predict?
a   = 384399000;
ecc = 0.0549;
r_horizons = a * (1 - ecc^2) / (1 + ecc * cos(nu_horizons));

% What r does your sim give?
M_sim  = deg2rad(135.49);
nu_sim = meanToTrue(M_sim, 0.0549);
r_sim  = a * (1 - ecc^2) / (1 + ecc * cos(nu_sim));

fprintf('r_horizons = %.0f km\n', r_horizons/1e3);
fprintf('r_sim      = %.0f km\n', r_sim/1e3);
fprintf('delta_r    = %.0f km\n', (r_sim - r_horizons)/1e3);