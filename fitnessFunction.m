function fitness = fitnessFunction(chromosome, params)
% -------------------------------------------------------------------------
% Function: fitnessFunction
% Purpose:
% Computes the fitness score for a given chromosome by evaluating the 
% trajectory's Delta-V and applying penalties for constraint violations 
% such as invalid sequences, missing planets, and excessive consecutive flybys.
%
% Inputs:
%   chromosome - Encoded trajectory [start_date, relative_dates, flyby_planets].
%   params     - Parameter structure containing mission-specific settings.
%
% Outputs:
%   fitness - Fitness score (lower is better). Infeasible trajectories are 
%             penalised with a fitness value of `Inf`.
%
% Dependencies:
%   computeTrajectory - Evaluates the trajectory's feasibility and Delta-V.
% -------------------------------------------------------------------------

    % Extract trajectory components from the chromosome
    start_date = chromosome(1);  % Start date (Julian)
    relativeDates = chromosome(2:params.maxFlybys+1);  % Relative flyby dates
    flybyPlanets = chromosome(params.maxFlybys+2:end); % Flyby planets (IDs)

    % Compute the trajectory
    [legs, totalDeltaV, isFeasible] = computeTrajectory(start_date, relativeDates, flybyPlanets, params);

    % Initialise fitness score
    if isFeasible
        % Base fitness is the total Delta-V of the trajectory
        fitness = totalDeltaV;

        % Penalty for excessive consecutive flybys of the same planet
        maxConsecutiveFlybys = 2; % Maximum allowed consecutive flybys
        consecutivePenalty = 0;
        for j = 1:length(flybyPlanets) - maxConsecutiveFlybys
            if all(flybyPlanets(j:j + maxConsecutiveFlybys - 1) == flybyPlanets(j))
                consecutivePenalty = consecutivePenalty + 300; % Penalty per occurrence
            end
        end

        % Penalty for missing required planets (Venus, Earth, Mars)
        requiredPlanets = [2, 3, 4]; % Venus (2), Earth (3), Mars (4)
        missingPlanets = setdiff(requiredPlanets, flybyPlanets);
        missingPlanetsPenalty = 300 * length(missingPlanets); % Penalty per missing planet

        % Penalty for invalid planet transition sequences
        invalidSequencePenalty = 0;
        for j = 1:length(flybyPlanets) - 1
            currentPlanet = flybyPlanets(j);
            nextPlanet = flybyPlanets(j + 1);
            % Disallow Venus-to-Mars or Mars-to-Venus transitions
            if (currentPlanet == 2 && nextPlanet == 4) || (currentPlanet == 4 && nextPlanet == 2)
                invalidSequencePenalty = invalidSequencePenalty + 300; % Penalty per invalid transition
            end
        end

        % Add penalties to the fitness score
        fitness = fitness + consecutivePenalty + missingPlanetsPenalty + invalidSequencePenalty;

    else
        % Assign infinite fitness for infeasible trajectories
        fitness = Inf;
    end
end
