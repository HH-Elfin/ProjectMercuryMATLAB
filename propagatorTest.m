clear; clc;

%% Constants
G = 6.674e-11; % m^3 kg^-1 s^-2

%% Time range (mirrors test_propagator.gd exactly)
t_start = 1747502700.0;                  % posixtime of 2026-05-17 18:25:00 UTC
t_end   = t_start + 365.25 * 86400.0;   % one Julian year later
steps   = 10;

%% Bodies
bodies = celestialBodies_propagatorTest();
bodies = linkParents(bodies, G);

%% Sample and print
for i = 1:steps
    t = t_start + (t_end - t_start) * (i - 1) / (steps - 1);

    bodies = updateBodyPositions(bodies, t);

    fprintf('--- t = %.0f ---\n', t);
    for j = 1:length(bodies)
        if bodies(j).name == "Sol"
            continue % Sol has no orbit — skip, matches Godot frame which holds only orbiting bodies
        end
        pos = bodies(j).position; % metres, heliocentric ecliptic
        fprintf('%s: (%.1f, %.1f, %.1f) m\n', bodies(j).name, pos(1), pos(2), pos(3));
    end
end