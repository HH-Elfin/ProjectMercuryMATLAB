function handles = plotSolarSystem(ax, bodies, handles, referenceBody, t)
    AU = 1.496e11;

    if nargin < 4 || isempty(referenceBody)
        referenceBody = 'Sol';
    end

    if nargin < 5
        t = [];
    end

    % Find the reference body's position (solar-frame)
    refPos = [0; 0; 0];
    for i = 1:length(bodies)
        if strcmp(bodies(i).name, referenceBody)
            if ~any(isnan(bodies(i).position))
                refPos = bodies(i).position;
            end
            break;
        end
    end

    for i = 1:length(bodies)
        if any(isnan(bodies(i).position))
            continue;
        end

        pos_AU = (bodies(i).position - refPos) / AU;
        name = bodies(i).name;

        % --- Orbit line ---
        if ~isempty(bodies(i).parent)
            parentPos = bodies(bodies(i).parent).position;
            orbitField = matlab.lang.makeValidName(name + "_orbit");
            baseField  = matlab.lang.makeValidName(name + "_orbit_base");

            if ~isfield(handles, orbitField)
                % First frame: compute and store base ellipse, create line
                basePts = computeOrbitEllipse(bodies(i));
                handles.(baseField) = basePts;
                displayPts = (basePts + parentPos - refPos) / AU;
                handles.(orbitField) = plot3(ax, ...
                    displayPts(1,:), displayPts(2,:), displayPts(3,:), ...
                    '-', 'Color', [bodies(i).colour, 0.3], ...
                    'LineWidth', 0.5, 'HandleVisibility', 'off');
            else
                % Subsequent frames: re-offset base points and update line
                basePts    = handles.(baseField);
                displayPts = (basePts + parentPos - refPos) / AU;
                handles.(orbitField).XData = displayPts(1,:);
                handles.(orbitField).YData = displayPts(2,:);
                handles.(orbitField).ZData = displayPts(3,:);
            end
        end

        % --- Body marker ---
        if isfield(handles, name)
            handles.(name).XData = pos_AU(1);
            handles.(name).YData = pos_AU(2);
            handles.(name).ZData = pos_AU(3);
            handles.(name).MarkerFaceColor = bodies(i).colour;
        else
            handles.(name) = plot3(ax, pos_AU(1), pos_AU(2), pos_AU(3), ...
                'o', 'MarkerSize', 6, 'DisplayName', name, ...
                'MarkerFaceColor', bodies(i).colour, ...
                'MarkerEdgeColor', 'none');
        end
    end

    if isfield(handles, 'timeTitle') && ~isempty(t)
        handles.timeTitle.String = formatSimTime(t);
    end

    if ~isfield(handles, 'initialised')
        grid(ax, 'on');
        axis(ax, 'equal');
        xlabel(ax, sprintf('X (AU) [origin: %s]', referenceBody));
        ylabel(ax, sprintf('Y (AU) [origin: %s]', referenceBody));
        zlabel(ax, sprintf('Z (AU) [origin: %s]', referenceBody));
        legend(ax, 'show');
        view(ax, 3);
        if ~isempty(t)
            handles.timeTitle = title(ax, formatSimTime(t), 'Color', [1 1 1]);
        end
        handles.initialised = true;
    end
end

function str = formatSimTime(t)
    str = datestr(datetime(t, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC'), ...
        'dd mmm yyyy HH:MM:SS UTC');
end