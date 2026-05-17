clear;clc;

%% Constants
G = 6.674e-11; % Gravitational constant, m^3 kg^-1 s^-2

%% Simulation parameters

% Points in time defined relative to the Unix epoch (1 Jan 1970, 00:00:00 UTC)
t_start = posixtime(datetime(2026, 05, 17, 18, 25, 0, 'TimeZone', 'UTC')); % simulation start date and time
t_end   = posixtime(datetime(2046, 05, 17, 18, 25, 0, 'TimeZone', 'UTC')); % simulation end date and time
timeWarpFactor = 1; % Simulation speed factor
realStepDuration = 1/30; % Real time between sim steps in seconds, essentially framerate
dt = timeWarpFactor * realStepDuration; % Each timestep is one second of real time, so at 1000x timewarp, dt is 1000 seconds
scaleAU = 0.1; % Radius of the display in astronomical units.
referenceBody = 'Earth'; % or 'Earth', 'Mars', etc.

%% Initialise

bodies = celestialBodies_sol();
bodies = linkParents(bodies, G);
t = t_start;
stepCount = 0;

%% Simulation

% Plot setup
fig = figure(1);  % Always use figure 1
clf(fig);         % Clear it completely, including all plot objects
ax = axes(fig);
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.ZLimMode = 'manual';
ax.Color = [0 0 0];
ax.GridColor = [1 1 1];
ax.XColor = [1 1 1];
ax.YColor = [1 1 1];
ax.ZColor = [1 1 1];
fig.Color = [0 0 0];
hold(ax, 'on');
plotHandles = struct();

uicontrol(fig, 'Style', 'pushbutton', 'String', 'Stop', ...
    'Position', [10 10 60 30], ...
    'Callback', @(~,~) setappdata(fig, 'stop', true));

setappdata(fig, 'stop', false);

while t < t_end && ishandle(fig)
    tic;

    % Log current step
    stepCount = stepCount + 1;
    fprintf('Step %d | t = %.0f | %s\n', stepCount, t, datetime(t, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC'));

    % Update celestial body positions to the correct one at the current timestep
    bodies = updateBodyPositions(bodies, t);
    
    % Visualisation
    if stepCount == 1
        cla(ax);
    end
    plotHandles = plotSolarSystem(ax, bodies, plotHandles, referenceBody, t);
    drawnow;

    % Fix limits
    if stepCount == 1
        axis(ax, [-scaleAU scaleAU -scaleAU scaleAU -scaleAU/4 scaleAU/4]);
        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';
        ax.ZLimMode = 'manual';
    end

    % Allow user to exit
    if getappdata(fig, 'stop')
        break;
    end

    % Wait for one second
    elapsed = toc;
    pauseTime = realStepDuration  - elapsed;
    if pauseTime > 0
        pause(pauseTime);
    end

    % Go to next timestep
    t = t + dt;
end