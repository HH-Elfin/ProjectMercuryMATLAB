function bodies = linkParents(bodies, G)
    for i = 1:length(bodies)
        if bodies(i).parentBody == ""
            bodies(i).parent = [];
            bodies(i).position = [0; 0; 0];
        else
            match = find([bodies.name] == bodies(i).parentBody);
            if isempty(match)
                error('Parent "%s" not found for body "%s"', bodies(i).parentBody, bodies(i).name);
            end
            bodies(i).parent = match;
        end

        % In linkParents, after resolving parent index:
        if ~isempty(bodies(i).parent)
            bodies(i).mu = G * (bodies(bodies(i).parent).mass + bodies(i).mass);
        end
        if bodies(i).name == "Luna"
            T_luna = 27.321582 * 86400;
            a_luna = (bodies(i).mu * (T_luna / (2*pi))^2)^(1/3);
            fprintf('Implied mean a = %.6e m\n', a_luna);
        end
    end
end