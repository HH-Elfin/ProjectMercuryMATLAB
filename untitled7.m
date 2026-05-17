clear;clc;
% What a do you get from the correct sidereal period?
G       = 6.674e-11;
M_earth = 5.97217e24;
M_luna  = 7.346e22;
mu      = G * (M_earth + M_luna);

T_sid   = 27.321582 * 86400;          % seconds
n_true  = 2*pi / T_sid;               % rad/s
a_true  = (mu / n_true^2)^(1/3);      % metres

fprintf('a_true = %.6e m  (%.3f km)\n', a_true, a_true/1e3);
fprintf('a_used = %.6e m  (%.3f km)\n', 3.796910736298334e8, 3.796910736298334e8/1e3);
fprintf('delta  = %.3f km\n', (a_true - 3.796910736298334e8)/1e3);