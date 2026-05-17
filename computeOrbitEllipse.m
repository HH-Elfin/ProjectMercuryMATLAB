function pts = computeOrbitEllipse(body, nPoints)
    if nargin < 2
        nPoints = 200;
    end

    a     = body.semiMajorAxis;
    e     = body.eccentricity;
    inc   = deg2rad(body.inclination);
    RAAN  = deg2rad(body.rightAscensionOfAscendingNode);
    omega = deg2rad(body.argumentOfPeriapsis);

    nu = linspace(0, 2*pi, nPoints);

    % Perifocal position
    r = a * (1 - e^2) ./ (1 + e * cos(nu));
    x = r .* cos(nu);
    y = r .* sin(nu);

    % Rotation matrix (same 3-1-3 as updateBodyPositions)
    R_omega = [cos(omega) -sin(omega) 0; sin(omega)  cos(omega) 0; 0 0 1];
    R_inc   = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
    R_RAAN  = [cos(RAAN) -sin(RAAN) 0; sin(RAAN)  cos(RAAN) 0; 0 0 1];
    R       = R_RAAN * R_inc * R_omega;

    % Apply rotation to all points at once
    pts = R * [x; y; zeros(1, nPoints)];
end