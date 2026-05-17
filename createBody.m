function body = createBody(options)
    arguments
        options.name string % Name of the body.
        options.mass (1,1) double {mustBePositive} % Mass of the body in kg.
        options.semiMajorAxis (1,1) double = NaN % Semi-major axis of the orbit in metres.
        options.eccentricity (1,1) double {mustBeInRange(options.eccentricity, 0, 1)} = 0 % Eccentricity of the orbit. Must be between 0 and 1.
        options.inclination (1,1) double = NaN % Inclination of the orbit in degrees.
        options.rightAscensionOfAscendingNode (1,1) double = NaN % Right ascension of ascending node of the orbit in degrees.
        options.argumentOfPeriapsis (1,1) double = NaN
        options.meanAnomaly (1,1) double = NaN % Mean anomaly of the body in degrees.
        options.isBarycentre (1,1) logical = false
        options.parentBody string = "" % Name of the parent body.
        options.position (3,1) double = [NaN; NaN; NaN] % Cartesian position (x,y,z). Calculated at runtime.
        options.velocity (3,1) double = [NaN; NaN; NaN] % Cartesian velocity (u,v,w). Calculated at runtime.
        options.mu (1,1) double = NaN % Gravitational parameter mu (G*M).
        options.colour (1,3) double {mustBeInRange(options.colour, 0, 1)} = [1, 1, 1] % Display colour for the planet.
        options.raanPrecessionRate (1,1) double = 0
        options.apsidalPrecessionRate (1,1) double = 0
    end

    body = options;
end