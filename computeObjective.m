function objective = computeObjective(start_date, flyby_dates, flyby_planets, params)
% -------------------------------------------------------------------------
% Function: computeObjective
% Purpose:
% Computes the objective value for trajectory optimisation by combining the 
% physical Delta-V and penalty contributions.
%
% Inputs:
%   start_date    - Start date of the mission (Julian date).
%   flyby_dates   - Relative dates of the flybys (days from start_date).
%   flyby_planets - Sequence of flyby planets (planet IDs).
%   params        - Parameter structure containing mission-specific settings.
%
% Outputs:
%   objective - Combined objective value (physical Delta-V + penalties).
% -------------------------------------------------------------------------

    % Compute physical Delta-V and penalty contributions
    [totalDeltaV, penalty] = objectiveFunction(start_date, flyby_dates, flyby_planets, params);

    % Combine physical Delta-V and penalty into a single objective value
    objective = totalDeltaV + penalty;

    % Debugging: Log the objective value
    fprintf('Objective = Physical Delta-V: %.3f km/s\n', objective);
end
