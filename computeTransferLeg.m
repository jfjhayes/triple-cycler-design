function [delta_v, leg_tof, success] = computeTransferLeg(departure_planet, arrival_planet, ...
    departure_date, arrival_date, params, nrev)

    % Default outputs for invalid solutions
    delta_v = Inf;
    success = false;

    % Compute time of flight (TOF)
    leg_tof = (arrival_date - departure_date) * 86400; % Convert days to seconds

    % Get state vectors for departure and arrival planets
    [~, r_departure, v_departure] = planet_elements_and_sv_jd(departure_planet, departure_date);
    [~, r_arrival, v_arrival] = planet_elements_and_sv_jd(arrival_planet, arrival_date);

    % Attempt short-period Lambert solution
    [~, ~, ~, ERROR, v_dep, v_arr] = lambertMR(r_departure', r_arrival', leg_tof, ...
        params.mu_sun, 0, nrev, 0, 0); % prograde only

    % Fallback: Try long-period solution if short-period fails
    if ERROR ~= 0
        [~, ~, ~, ERROR, v_dep, v_arr] = lambertMR(r_departure', r_arrival', leg_tof, ...
            params.mu_sun, 0, nrev, 1, 0); % Long-period mode
    end

    % If both solutions fail, return invalid
    if ERROR ~= 0
        return; % neat error handling
    end

    % Compute delta-v for the transfer
    v_inf_departure = norm(v_dep - v_departure); % Hyperbolic excess at departure
    v_inf_arrival = norm(v_arr - v_arrival);     % Hyperbolic excess at arrival
    delta_v = v_inf_departure + v_inf_arrival;

    % Mark solution as valid
    success = true;
end
