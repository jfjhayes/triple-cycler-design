function [legs, totalDeltaV, isFeasible, logData] = computeTrajectory(start_date, relativeDates, flybyPlanets, params)
% -------------------------------------------------------------------------
% Function: computeTrajectory
% Purpose:
% Computes the trajectory for a given start date, relative flyby dates, 
% and flyby planets. It evaluates each transfer leg, computes Delta-V, and 
% logs detailed trajectory data.
%
% Inputs:
%   start_date    - Start date of the mission (Julian date).
%   relativeDates - Relative flyby dates (days from start_date).
%   flybyPlanets  - Sequence of flyby planets (planet IDs).
%   params        - Parameter structure containing mission-specific settings.
%
% Outputs:
%   legs         - Array of transfer legs and their properties.
%   totalDeltaV  - Total Delta-V required for the trajectory (km/s).
%   isFeasible   - Boolean indicating if the trajectory is feasible.
%   logData      - Struct array containing detailed trajectory logs.
%
% Dependencies:
%   computeTransferLeg    - Evaluates a single transfer leg.
%   planet_elements_and_sv_jd - Computes planetary state vectors.
%   computeLegGAM         - Checks gravity assist feasibility.
% -------------------------------------------------------------------------

    % Initialise outputs
    legs = [];
    totalDeltaV = 0;
    isFeasible = true;

    % Initialise logData structure
    logData = struct('leg', [], 'dep_date_jd', [], 'arr_date_jd', [], ...
                     'dep_date_cal', [], 'arr_date_cal', [], ...
                     'dep_planet', [], 'arr_planet', [], ...
                     'v_inf_in', [], 'v_inf_out', [], ...
                     'delta_v', [], 'tof', [], ...
                     'dep_r', [], 'dep_v', [], ...
                     'arr_r', [], 'arr_v', []);

    % Add dummy start and end planets
    dummy_planet = params.startPlanet; % Start and end planet (default: Earth)
    planets = [dummy_planet, flybyPlanets(:)', dummy_planet];
    dates = [start_date, start_date + cumsum(relativeDates), start_date + 2 * params.synodicPeriod];

    % Ensure dates are strictly increasing
    for i = 2:length(dates)
        if dates(i) <= dates(i - 1)
            dates(i) = dates(i - 1) + 1e-6; % Adjust to prevent overlapping dates
        end
    end

    % Loop through each transfer leg
    for i = 1:(length(planets) - 1)
        dep_planet = planets(i);       % Departure planet
        arr_planet = planets(i + 1);  % Arrival planet
        dep_date = dates(i);           % Departure date (Julian)
        arr_date = dates(i + 1);       % Arrival date (Julian)

        % Compute the transfer leg
        disp(['Computing leg ', num2str(i), ' of ', num2str(length(planets) - 1)]);
        [delta_v, leg_tof, success, v_departure, v_arrival] = ...
            computeTransferLeg(dep_planet, arr_planet, dep_date, arr_date, params);

        % Log data for the current leg
        logData(i).leg = i;
        logData(i).dep_date_jd = dep_date;
        logData(i).arr_date_jd = arr_date;
        logData(i).dep_date_cal = datetime(dep_date, 'ConvertFrom', 'juliandate', 'Format', 'yyyy-MM-dd HH:mm:ss');
        logData(i).arr_date_cal = datetime(arr_date, 'ConvertFrom', 'juliandate', 'Format', 'yyyy-MM-dd HH:mm:ss');
        logData(i).dep_planet = dep_planet;
        logData(i).arr_planet = arr_planet;
        logData(i).v_inf_in = norm(v_departure);
        logData(i).v_inf_out = norm(v_arrival);
        logData(i).delta_v = delta_v;
        logData(i).tof = leg_tof / 86400; % Convert TOF to days

        % Log state vectors for departure and arrival
        [~, dep_r, dep_v] = planet_elements_and_sv_jd(dep_planet, dep_date);
        [~, arr_r, arr_v] = planet_elements_and_sv_jd(arr_planet, arr_date);
        logData(i).dep_r = dep_r;
        logData(i).dep_v = dep_v;
        logData(i).arr_r = arr_r;
        logData(i).arr_v = arr_v;

        % Check feasibility of the transfer leg
        if ~success
            isFeasible = false;
            totalDeltaV = Inf; % Mark as infeasible
            return;
        end

        % Accumulate Delta-V for the leg
        totalDeltaV = totalDeltaV + delta_v;
    end

    % Perform gravity assist feasibility checks
    [isFeasible, assistDeltaV] = computeLegGAM(legs, planets, params);
    if isFeasible
        totalDeltaV = totalDeltaV + assistDeltaV;
    else
        totalDeltaV = Inf; % Mark as infeasible if gravity assist is invalid
    end
end
