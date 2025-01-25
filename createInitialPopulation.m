function population = createInitialPopulation(popSize, params)
% -------------------------------------------------------------------------
% Function: createInitialPopulation
% Purpose:
% Generates the initial population for the Genetic Algorithm (GA) by 
% randomly creating chromosomes that satisfy the mission constraints.
%
% Inputs:
%   popSize - Number of individuals in the population.
%   params  - Parameter structure containing mission-specific settings:
%             - start_date: Earliest possible mission start date (Julian).
%             - end_date: Latest possible mission end date (Julian).
%             - maxFlybys: Maximum number of flybys allowed.
%
% Outputs:
%   population - Array of structures representing the population:
%                - chromosome: Encoded trajectory [start_date, relative_dates, flyby_planets].
%                - start_date: Randomly selected start date (Julian).
%                - flyby_dates_relative: Relative flyby dates (days from start_date).
%                - flyby_planets: Sequence of flyby planets (IDs).
%
% Dependencies:
%   None.
% -------------------------------------------------------------------------

    % Constants
    synodicPeriod = 4675.2; % Mission duration (12.8 years) in days
    maxFlybys = params.maxFlybys; % Maximum number of flybys allowed

    % Preallocate the population structure
    population = repmat(struct('chromosome', [], 'start_date', [], ...
                               'flyby_dates_relative', [], 'flyby_planets', []), popSize, 1);

    for i = 1:popSize
        % Generate a random start date within the allowable range
        t_start = params.start_date + floor(rand * (params.end_date - params.start_date));

        % Define the mission end date
        t_end = t_start + synodicPeriod;

        % Generate random relative flyby dates
        valid = false;
        while ~valid
            % Generate sorted random intervals
            flyby_intervals = sort(rand(1, maxFlybys)); 
            flyby_intervals = flyby_intervals / sum(flyby_intervals) * synodicPeriod; % Scale to synodic period
            % Ensure the total duration is valid
            if sum(flyby_intervals) < synodicPeriod
                valid = true;
            end
        end

        % Compute cumulative flyby times
        t_flybys = cumsum(flyby_intervals);

        % Compute relative dates as intervals
        relativeDates = diff([0, t_flybys]);

        % Adjust relative dates if the sum exceeds the synodic period
        if sum(relativeDates) >= synodicPeriod
            warning('Rescaling relative dates to fit within mission duration.');
            relativeDates = relativeDates / sum(relativeDates) * (synodicPeriod - rand * 100);
        end

        % Generate random flyby planets (Venus, Earth, Mars)
        flybyPlanets = randi([2, 4], 1, maxFlybys);

        % Construct the chromosome
        chromosome = [t_start, relativeDates, flybyPlanets];

        % Populate the population structure
        population(i).chromosome = chromosome;
        population(i).start_date = t_start;
        population(i).flyby_dates_relative = relativeDates;
        population(i).flyby_planets = flybyPlanets;

        % Debugging: Display initial population details
        % disp('Generated Individual:');
        % disp(population(i));
    end
end
