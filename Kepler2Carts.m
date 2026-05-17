function [X, Y, Z, Vx, Vy, Vz] = Kepler2Carts(a, ecc, inc, w, nu, RAAN, mu)
%-------------------------------------------------------------------------%
%
%   USAGE: Conversion of Keplerian orbital elements to Cartesian position
%          and velocity in the parent body's inertial reference frame.
%
%   AUTHOR: Adapted from Thameur Chebbi (PhD), chebbythamer@gmail.com
%           Original: 01 Oct 2020
%           Modified: generalised for arbitrary parent bodies
%
%   DESCRIPTION:
%       Converts classical orbital elements to Cartesian state vectors
%       (position and velocity) in the perifocal frame, then rotates into
%       the parent body's inertial frame via the standard 3-1-3 Euler
%       rotation sequence (RAAN, inc, w).
%
%       Units are caller-defined. mu, a, and outputs will be consistent
%       with whatever unit system is passed in. Recommended conventions:
%           Solar system scale  — AU, AU/day,  mu in AU^3/day^2
%           Planetary scale     — km, km/s,    mu in km^3/s^2
%
%   INPUT:
%       a       Semi-major axis              (caller's length unit)
%       ecc     Eccentricity                 (dimensionless)
%       inc     Inclination                  (rad)
%       w       Argument of periapsis        (rad)
%       nu      True anomaly                 (rad)
%       RAAN    Right ascension of asc. node (rad)
%       mu      Gravitational parameter of   (length^3 / time^2)
%               parent body (G * M_parent)
%
%   OUTPUT:
%       [X, Y, Z]       Position vector      (caller's length unit)
%       [Vx, Vy, Vz]    Velocity vector      (caller's length unit / time)
%
%   NOTE:
%       nu (true anomaly) must be computed by the caller from the mean
%       anomaly M(t) via Kepler's equation. This function does not solve
%       Kepler's equation internally.
%
%-------------------------------------------------------------------------%

%% Semi-latus rectum and radial distance
p   = a * (1 - ecc^2);
r_0 = p / (1 + ecc * cos(nu));

%% Perifocal frame — position
x = r_0 * cos(nu);
y = r_0 * sin(nu);

%% Perifocal frame — velocity
Vx_ = -(mu / p)^0.5 * sin(nu);
Vy_ =  (mu / p)^0.5 * (ecc + cos(nu));

%% Rotation matrix coefficients (3-1-3: RAAN -> inc -> w)
cO = cos(RAAN);  sO = sin(RAAN);
ci = cos(inc);   si = sin(inc);
cw = cos(w);     sw = sin(w);

R11 =  cO*cw - sO*sw*ci;
R12 = -cO*sw - sO*cw*ci;
R21 =  sO*cw + cO*sw*ci;
R22 = -sO*sw + cO*cw*ci;
R31 =  sw*si;
R32 =  cw*si;

%% Rotate position into inertial frame
X = R11*x + R12*y;
Y = R21*x + R22*y;
Z = R31*x + R32*y;

%% Rotate velocity into inertial frame
Vx = R11*Vx_ + R12*Vy_;
Vy = R21*Vx_ + R22*Vy_;
Vz = R31*Vx_ + R32*Vy_;

end