function [leg_dv, leg_tof, leg_data] = propagateSegment(dates, b_mag, b_angle, params)

    % Debug output
    fprintf('Processing dates: %.2f to %.2f\n', dates(1), dates(2));
    fprintf('Raw dates: %.10f, %.10f\n', dates(1), dates(2));

    % Convert Julian dates
    [y1,m1,d1] = julian2calendar(dates(1));
    [y2,m2,d2] = julian2calendar(dates(2));
    
    fprintf('Date 1: %d-%d-%d\n', y1, m1, d1);
    fprintf('Date 2: %d-%d-%d\n', y2, m2, d2);
    
    % Get planet states with debug output
    [elements1, r1, v1, jd1] = planet_elements_and_sv(params.planet1, y1, m1, d1, 0, 0, 0);
    [elements2, r2, v2, jd2] = planet_elements_and_sv(params.planet2, y2, m2, d2, 0, 0, 0);
    
    fprintf('Planet 1 state (ID %d):\nr: [%.2e %.2e %.2e]\nv: [%.2e %.2e %.2e]\n', ...
        params.planet1, r1(1), r1(2), r1(3), v1(1), v1(2), v1(3));
    fprintf('Planet 2 state (ID %d):\nr: [%.2e %.2e %.2e]\nv: [%.2e %.2e %.2e]\n', ...
        params.planet2, r2(1), r2(2), r2(3), v2(1), v2(2), v2(3));
    
    % Ensure correct dimensions 
    r1 = reshape(r1, [3,1]);
    r2 = reshape(r2, [3,1]);
    v1 = reshape(v1, [3,1]);
    v2 = reshape(v2, [3,1]);
    
    leg_tof = (jd2 - jd1) * 86400;

    % Set multi-rev parameters
    nrev = 0;
    orbit_type = 0;
    if isfield(params, 'nrev')
        nrev = params.nrev;
    end
    if isfield(params, 'orbit_type')
        orbit_type = params.orbit_type;
    end

    % Validate state vectors
    if any(isnan([r1; r2; v1; v2])) || leg_tof <= 0
        error('Invalid state vectors or transfer time');
    end
    if norm(r1) < 1e7 || norm(r2) < 1e7  % km
        error('Position vectors too close to sun');
    end

    % Lambert solution with multi-rev options
    [A, P, E, ERROR, v_dep, v_arr] = lambertMR(r1', r2', leg_tof, ...
    params.mu_sun, orbit_type, nrev, 0, 2); % Enable full warnings

    % If short-period fails, try long-period
    if ERROR ~= 0
        [A, P, E, ERROR, v_dep, v_arr] = lambertMR(r1', r2', leg_tof, ...
        params.mu_sun, orbit_type, nrev, 1, 2);
    end

    if ERROR ~= 0
        error('Lambert solver failed with error code %d', ERROR);
        leg_dv = Inf;
        leg_tof = Inf;
        leg_data = struct('altitude', Inf, 'v_in', zeros(3,1), 'v_out', zeros(3,1), ...
                         'v_planet', v1, 'r_planet', r1, 'r_target', r2, 'v_target', v2);
        return
    end
    
    fprintf('Lambert SMA %g\n', A);
    fprintf('Lambert SLR %g\n', P);
    fprintf('Lambert E %g\n', E);

    leg_data = struct('altitude', 0, 'turn_angle', 0, ...
                     'v_in', zeros(3,1), 'v_out', zeros(3,1), ...
                     'v_planet', v1, 'r_planet', r1, ...
                     'r_target', r2, 'v_target', v2, ...
                     'feasible', false);

    planet_data = planet_properties(params.planet1);
    [v_out, alt, turn] = computeGravityAssist(v_dep', v1, planet_data, b_mag, b_angle);

    % Store results
    leg_data.altitude = alt;
    leg_data.turn_angle = turn;
    leg_data.v_in = v_dep';
    leg_data.v_out = v_out;
    leg_data.feasible = true;
    
    % Validate and store trajectory info
    validateAndStoreTrajInfo(leg_data, v_dep', v_arr', r1, r2, v1, v2, params.mu_sun);
    
    % Calculate deltaV based on feasibility
    if ~leg_data.feasible
        v_inf_delta = abs(norm(leg_data.v_out - leg_data.v_planet) - ...
                         norm(leg_data.v_in - leg_data.v_planet));
        leg_dv = v_inf_delta;
    else
        leg_dv = 0;
    end
end