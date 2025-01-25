function [totalDeltaV, penalty] = objectiveFunction(start_date, flyby_dates, flyby_planets, params)
% -------------------------------------------------------------------------
% Function: objectiveFunction
% Purpose:
% Evaluates the physical Delta-V and penalty for a trajectory over two 
% passes (for improved accuracy).
%
% Inputs:
%   start_date    - Start date of the mission (Julian date).
%   flyby_dates   - Relative dates of the flybys (days from start_date).
%   flyby_planets - Sequence of flyby planets (planet IDs).
%   params        - Parameter structure containing mission-specific settings.
%
% Outputs:
%   totalDeltaV - Total physical Delta-V of the trajectory (km/s).
%   penalty     - Total penalty due to constraint violations.
% -------------------------------------------------------------------------

    % First pass: Compute Delta-V, penalty, and final velocity
    fprintf('--- First Pass ---\n');
    [totalDeltaV, penalty, final_velocity] = analyzeTrajectory(start_date, flyby_dates, flyby_planets, params, []);

    % Second pass: Refine Delta-V and penalty using final velocity from the first pass
    fprintf('--- Second Pass ---\n');
    [totalDeltaV, penalty] = analyzeTrajectory(start_date, flyby_dates, flyby_planets, params, final_velocity);
end


function [totalDeltaV, penalty, final_velocity] = analyzeTrajectory(start_date, flyby_dates, flyby_planets, params, initial_velocity)
% -------------------------------------------------------------------------
% Function: analyzeTrajectory
% Purpose:
% Computes the total Delta-V and penalties for a given trajectory by 
% analysing each leg of the transfer, considering Lambert transfers and 
% gravity assist constraints.
%
% Inputs:
%   start_date      - Start date of the mission (Julian date).
%   flyby_dates     - Relative dates of the flybys (days from start_date).
%   flyby_planets   - Sequence of flyby planets (planet IDs).
%   params          - Parameter structure containing mission-specific settings.
%   initial_velocity - Initial velocity vector for the first transfer leg (optional).
%
% Outputs:
%   totalDeltaV   - Total physical Delta-V of the trajectory (km/s).
%   penalty       - Total penalty due to constraint violations.
%   final_velocity - Final velocity vector at the end of the trajectory.
% -------------------------------------------------------------------------

    % Define start and end planets
    start_planet = 4;  % Mars (default for outbound)
    end_planet = 4;    % Mars (wrap-around trajectory)

    % Compute absolute dates from relative flyby dates
    absolute_dates = start_date + cumsum([0, flyby_dates, params.end_date - start_date]);

    % Extend flyby_planets to include start and end planets
    extended_planets = [start_planet, flyby_planets, end_planet];
    disp('Extended Planet Sequence:');
    disp(extended_planets);

    % Initialise outputs
    totalDeltaV = 0;
    penalty = 0;
    final_velocity = [];

    % Iterate through all trajectory legs
    for i = 1:length(extended_planets) - 1
        current_planet = extended_planets(i);
        next_planet = extended_planets(i + 1);

        % Get planet state vectors
        [~, r1, v1, mu1] = planet_elements_and_sv_jd(current_planet, absolute_dates(i));
        [~, r2, v2, ~] = planet_elements_and_sv_jd(next_planet, absolute_dates(i + 1));

        % Use initial velocity if provided for the first leg
        if i == 1 && ~isempty(initial_velocity)
            v1 = initial_velocity;
        end

        % Solve Lambert problem
        tof_seconds = (absolute_dates(i + 1) - absolute_dates(i)) * 86400;  % Time of flight in seconds
        [~, ~, ~, error, vi, vf] = lambertMR(r1, r2, tof_seconds, params.mu_sun, 0, 0, 0, 0);
        if error ~= 0
            fprintf('Lambert solution failed for leg %d-%d\n', current_planet, next_planet);
            return;
        end

        % Compute Delta-V for the leg
        deltaV_leg = norm(vi - v1) + norm(vf - v2);
        totalDeltaV = totalDeltaV + deltaV_leg;

        % Gravity assist constraints and penalties
        if i < length(extended_planets) - 1
            v_inf_in = vi - v1;  % Incoming hyperbolic excess velocity
            v_inf_out = vf - v2; % Outgoing hyperbolic excess velocity

            % Compute periapsis radius
            periapsis_radius = computePeriapsisRadius(v_inf_in, v_inf_out, mu1);

            % Apply penalties for invalid periapsis
            planet_radius = params.planets(current_planet - 1).radius;
            min_periapsis = planet_radius + params.min_altitude;
            max_periapsis = planet_radius + 100000;  % Maximum allowable periapsis radius
            if periapsis_radius < min_periapsis || periapsis_radius > max_periapsis
                penalty = penalty + abs(periapsis_radius - min_periapsis) * 1e3;
            end
        end

        % Save final velocity for wrap-around trajectory
        if i == length(flyby_planets) - 1
            final_velocity = vf;
        end
    end

    % Debugging: Log results
    fprintf('Total Physical Delta-V = %.3f km/s\n', totalDeltaV);
end


