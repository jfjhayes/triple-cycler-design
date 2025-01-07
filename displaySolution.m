function displaySolution(solution)
    fprintf('\nCycler Solution Summary:\n');
    fprintf('------------------------\n');
    
    for i = 1:solution.num_segments
        [y,m,d] = julian2calendar(solution.dates(i));
        fprintf('\nSegment %d:\n', i);
        fprintf('Departure: %s %d, %d from %s\n', ...
            month_name(m), d, y, planet_name(getPlanetID(i)));
            
        [y,m,d] = julian2calendar(solution.dates(i+1));
        fprintf('Arrival: %s %d, %d at %s\n', ...
            month_name(m), d, y, planet_name(getPlanetID(i+1)));
        
        fprintf('TOF: %.1f days\n', solution.tof(i));
        fprintf('Flyby alt: %.1f km\n', solution.altitudes(i));
        fprintf('V∞: %.2f km/s\n', solution.vinf(i));
    end
end

function solution = processSolution(x_best, constants, constraints)
    [dates, b_params] = decodeChromosome(x_best);
    num_segments = size(dates,1)-1;
    
    solution = struct(...
        'dates', dates,...
        'b_params', b_params,...
        'total_dv', 0,...
        'num_segments', num_segments,...
        'tof', zeros(num_segments,1),...
        'altitudes', zeros(num_segments,1),...
        'vinf', zeros(num_segments,1),...
        'planets', cell(num_segments+1,1));
        
    for i = 1:num_segments
        p_id = getPlanetID(i);
        solution.planets{i} = planet_name(p_id);
        
        seg_params = struct(...
            'planet1', p_id,...
            'planet2', getPlanetID(i+1),...
            'mu_sun', constants.mu_sun);
            
        [dv, tof, leg_data] = propagateSegment(...
            dates(i:i+1,:), b_params(i,1), b_params(i,2), seg_params);
            
        solution.total_dv = solution.total_dv + dv;
        solution.tof(i) = tof/86400;  % Convert to days
        solution.altitudes(i) = leg_data.altitude;
        solution.vinf(i) = norm(leg_data.v_in - leg_data.v_planet);
    end
    
    solution.planets{end} = planet_name(getPlanetID(num_segments+1));
end

function p_id = getPlanetID(segment_num)
    sequence = [3,2,3,2,3,4];  % Earth-Venus-Earth-Venus-Earth-Mars
    p_id = sequence(mod(segment_num-1,length(sequence))+1);
end