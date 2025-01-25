function [isFeasible, totalDeltaV] = computeLegGAM(legs, planets, params)
% -------------------------------------------------------------------------
% Function: computeLegGAM
% Purpose:
% Evaluates the feasibility of gravity assists for a trajectory by computing 
% hyperbolic excess velocities and turning angles. Also calculates Delta-V 
% for Earth/Mars circular orbit transfers if applicable.
%
% Inputs:
%   legs    - Array of trajectory legs containing arrival and departure velocities.
%   planets - Sequence of planet IDs for the trajectory.
%   params  - Parameter structure containing mission-specific settings:
%             - min_altitude: Minimum allowable altitude for gravity assists.
%
% Outputs:
%   isFeasible   - Boolean indicating whether all gravity assists are feasible.
%   totalDeltaV  - Total Delta-V for gravity assists and circular orbit transfers.
%
% Dependencies:
%   planet_properties - Retrieves physical properties of planets (e.g., radius, mu).
% -------------------------------------------------------------------------

    % Initialise outputs
    isFeasible = true;
    totalDeltaV = 0;
    numLegs = length(legs);

    % Loop through trajectory legs, skipping dummy flybys
    for i = 1:numLegs
        if planets(i) == planets(end) % Skip the dummy flyby (final planet)
            continue;
        end

        % Retrieve planet properties
        planet_id = planets(i);
        planet = planet_properties(planet_id);

        % Compute hyperbolic excess velocities
        v_inf_in = legs(i).v_arr - planet.mu / planet.radius; % Incoming excess velocity
        v_inf_out = legs(i).v_dep - planet.mu / planet.radius; % Outgoing excess velocity

        % Calculate turning angle (delta) and maximum allowable turning angle (delta_max)
        delta = acos(dot(v_inf_in, v_inf_out) / (norm(v_inf_in) * norm(v_inf_out))); % Turning angle
        delta_max = 2 * asin(planet.mu / (planet.mu + params.min_altitude * norm(v_inf_in)^2)); % Max turning angle

        % Check gravity assist feasibility
        if delta > delta_max
            disp('Gravity assist infeasible: Turning angle exceeds maximum.');
            isFeasible = false;
            return; % Exit early if gravity assist is not feasible
        end

        % Compute taxi Delta-V for Earth/Mars circular orbit transfers
        if planet_id == 3 || planet_id == 4 % Earth (3) or Mars (4)
            legs(i).v_inf_avg = (norm(v_inf_in) + norm(v_inf_out)) / 2; % Average excess velocity
            legs(i).taxi_delta_v_in = norm(v_inf_in) / sqrt(2); % Incoming circularisation Delta-V
            legs(i).taxi_delta_v_out = norm(v_inf_out) / sqrt(2); % Outgoing escape Delta-V
            totalDeltaV = totalDeltaV + legs(i).taxi_delta_v_in + legs(i).taxi_delta_v_out; % Accumulate Delta-V
        end
    end
end
