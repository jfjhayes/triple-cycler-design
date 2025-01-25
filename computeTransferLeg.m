function [delta_v, leg_tof, success, v_departure, v_arrival] = computeTransferLeg(departure_planet, arrival_planet, ...
    departure_date, arrival_date, params)
% -------------------------------------------------------------------------
% Function: computeTransferLeg
% Purpose:
% Computes the Delta-V and time of flight (TOF) for a single transfer leg 
% between two planets using Lambert's problem solver. If the short-period 
% solution fails, the function attempts a fallback with the long-period mode.
%
% Inputs:
%   departure_planet - Planet ID for the departure planet.
%   arrival_planet   - Planet ID for the arrival planet.
%   departure_date   - Departure date (Julian).
%   arrival_date     - Arrival date (Julian).
%   params           - Parameter structure containing mission-specific settings.
%
% Outputs:
%   delta_v    - Total Delta-V for the transfer leg (km/s).
%   leg_tof    - Time of flight for the transfer leg (seconds).
%   success    - Boolean indicating whether the transfer leg is valid.
%   v_departure - Velocity vector at departure (km/s).
%   v_arrival   - Velocity vector at arrival (km/s).
%
% Dependencies:
%   planet_elements_and_sv_jd - Computes state vectors for planets.
%   lambertMR                 - Solves Lambert's problem for orbital transfer.
% -------------------------------------------------------------------------

    % Initialise outputs for invalid solutions
    delta_v = Inf;
    success = false;

    % Compute time of flight (TOF) in seconds
    leg_tof = (arrival_date - departure_date) * 86400;

    % Retrieve state vectors for departure and arrival planets
    [~, r_departure, v_departure] = planet_elements_and_sv_jd(departure_planet, departure_date);
    [~, r_arrival, v_arrival] = planet_elements_and_sv_jd(arrival_planet, arrival_date);

    % Debugging: Display state vectors for the Lambert solver
    disp('State Vectors for Lambert Solver:');
    disp(['r_departure: ', num2str(r_departure)]);
    disp(['v_departure: ', num2str(v_departure)]);
    disp(['r_arrival: ', num2str(r_arrival)]);
    disp(['v_arrival: ', num2str(v_arrival)]);

    % Attempt short-period Lambert solution
    [~, ~, ~, ERROR, v_dep, v_arr] = lambertMR(r_departure', r_arrival', leg_tof, ...
        params.mu_sun, 0, 0, 0, 0); % Prograde, short-period solution

    % Fallback: Try long-period Lambert solution if short-period fails
    if ERROR ~= 0
        [~, ~, ~, ERROR, v_dep, v_arr] = lambertMR(r_departure', r_arrival', leg_tof, ...
            params.mu_sun, 0, 0, 1, 0); % Prograde, long-period solution
    end

    % If both solutions fail, return invalid
    if ERROR ~= 0
        disp('Lambert Solver Failed:');
        disp(['ERROR Code: ', num2str(ERROR)]);
        return; % Terminate with invalid state
    end

    % Debugging: Display Lambert solver outputs
    disp('Lambert Solver Outputs:');
    disp(['v_dep (Spacecraft): ', num2str(v_dep)]);
    disp(['v_arr (Spacecraft): ', num2str(v_arr)]);

    % Compute Delta-V for the transfer
    v_inf_departure = norm(v_dep - v_departure); % Hyperbolic excess velocity at departure
    v_inf_arrival = norm(v_arr - v_arrival);     % Hyperbolic excess velocity at arrival
    delta_v = v_inf_departure + v_inf_arrival;

    % Mark the transfer leg as valid
    success = true;
end
