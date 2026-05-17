function handles = plotSolarSystem(ax, bodies, handles, referenceBody)
    AU = 1.496e11;

    % Default reference body is Sol
    if nargin < 4 || isempty(referenceBody)
        referenceBody = 'Sol';
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

        % Express position relative to the reference body
        pos_AU = (bodies(i).position - refPos) / AU;
        name = bodies(i).name;

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

    if ~isfield(handles, 'initialised')
        grid(ax, 'on');
        axis(ax, 'equal');
        xlabel(ax, sprintf('X (AU) [origin: %s]', referenceBody));
        ylabel(ax, sprintf('Y (AU) [origin: %s]', referenceBody));
        zlabel(ax, sprintf('Z (AU) [origin: %s]', referenceBody));
        legend(ax, 'show');
        view(ax, 3);
        handles.initialised = true;
    end
end