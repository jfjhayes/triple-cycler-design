function solution = processSolution(x_best, constants, constraints)
    [dates, b_params] = decodeChromosome(x_best);
    num_segments = length(dates) - 1;
    
    solution = struct(...
        'dates', dates, ...
        'b_params', b_params, ...
        'total_dv', 0, ...
        'num_segments', num_segments, ...
        'tof', zeros(num_segments,1), ...
        'altitudes', zeros(num_segments,1), ...
        'vinf', zeros(num_segments,1));
    
    for i = 1:num_segments
        seg_params = struct(...
            'planet1', getPlanetID(i), ...
            'planet2', getPlanetID(i+1), ...
            'mu_sun', constants.mu_sun);
            
        [dv, tof, leg_data] = propagateSegment(...
            [dates(i), dates(i+1)], b_params(i,1), b_params(i,2), seg_params);
            
        solution.total_dv = solution.total_dv + dv;
        solution.tof(i) = tof/86400;  % Convert to days
        solution.altitudes(i) = leg_data.altitude;
        solution.vinf(i) = norm(leg_data.v_in - leg_data.v_planet);
    end
end