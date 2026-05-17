function system = updateBodyPositions(bodies, t)
    system = bodies;

    % Build a processing order that guarantees parents are resolved before
    % children. parent is a numeric index (or [] for roots) set by linkParents.
    order = resolveOrder(system);

    for idx = 1:length(order)
        i = order(idx);

        if isempty(system(i).parent)
            continue; % Root body (e.g. Sol) — position already set by linkParents
        end

        % Orbital elements
        a     = system(i).semiMajorAxis;
        e     = system(i).eccentricity;
        inc   = deg2rad(system(i).inclination);
        RAAN  = deg2rad(system(i).rightAscensionOfAscendingNode);
        omega = deg2rad(system(i).argumentOfPeriapsis);
        M0    = deg2rad(system(i).meanAnomaly);
        mu    = system(i).mu;
        
        % Apply apsidal and nodal precession if specified (e.g. Luna)
        % Rates are deg/day; t is seconds since Unix epoch
        RAAN  = RAAN  + deg2rad(system(i).raanPrecessionRate  / 86400) * t;
        omega = omega + deg2rad(system(i).apsidalPrecessionRate / 86400) * t;

        % Step 1: Propagate mean anomaly
        n = sqrt(mu / a^3);         % Mean motion (rad/s)
        M = mod(M0 + n * t, 2*pi);  % Mean anomaly at time t

        % Step 2: Solve Kepler's equation for eccentric anomaly (Newton-Raphson)
        E = M;
        for iter = 1:100
            dE = (M - E + e * sin(E)) / (1 - e * cos(E));
            E  = E + dE;
            if abs(dE) < 1e-10
                break;
            end
        end

        % Step 3: True anomaly
        nu = 2 * atan2(sqrt(1+e) * sin(E/2), sqrt(1-e) * cos(E/2));

        % Step 4: Position in perifocal frame
        r      = a * (1 - e * cos(E));
        r_peri = r * [cos(nu); sin(nu); 0];

        % Step 5: Rotate to ecliptic frame
        R_omega = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
        R_inc   = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
        R_RAAN  = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
        R       = R_RAAN * R_inc * R_omega;

        r_ecliptic = R * r_peri;

        % Step 6: Add parent's ecliptic position.
        % parent is a numeric index into system, guaranteed resolved already.
        parentIdx = system(i).parent;
        system(i).position = r_ecliptic + system(parentIdx).position;
    end
end

% -------------------------------------------------------------------------
function order = resolveOrder(system)
    % Topological sort: returns body indices in an order where every parent
    % is processed before its children.
    % parent is [] for roots, or a numeric index set by linkParents.
    n      = length(system);
    order  = zeros(1, n);
    placed = false(1, n);
    count  = 0;

    % Place roots first
    for i = 1:n
        if isempty(system(i).parent)
            count        = count + 1;
            order(count) = i;
            placed(i)    = true;
        end
    end

    % Iteratively place bodies whose parent is already placed
    maxIter = n * n;
    iter    = 0;
    while count < n && iter < maxIter
        iter = iter + 1;
        for i = 1:n
            if placed(i)
                continue;
            end
            parentIdx = system(i).parent;
            if placed(parentIdx)
                count        = count + 1;
                order(count) = i;
                placed(i)    = true;
            end
        end
    end

    if count < n
        warning('updateBodyPositions:unresolvedParents', ...
            'Could not resolve parent order for all bodies. Check parentBody fields.');
        for i = 1:n
            if ~placed(i)
                count        = count + 1;
                order(count) = i;
            end
        end
    end
end