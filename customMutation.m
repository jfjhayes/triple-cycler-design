function mutatedChromosome = customMutation(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation, params)
% -------------------------------------------------------------------------
% Function: customMutation
% Purpose:
% Custom mutation function for the Genetic Algorithm (GA) to introduce 
% variations in chromosomes while enforcing constraints such as mission 
% duration and valid flyby planets.
%
% Inputs:
%   parents        - Indices of the parent chromosomes selected for mutation.
%   options        - Options structure for the GA.
%   nvars          - Number of variables in each chromosome.
%   FitnessFcn     - Fitness function for the GA.
%   state          - State structure of the GA at the current generation.
%   thisScore      - Scores of the population at the current generation.
%   thisPopulation - Current population of chromosomes.
%   params         - Parameter structure containing mission-specific settings.
%
% Outputs:
%   mutatedChromosome - Mutated version of the parent chromosome.
%
% Dependencies:
%   None.
% -------------------------------------------------------------------------

    % Mission duration (in days)
    synodicPeriod = params.synodicPeriod; % 12.8 years

    % Ensure parents are selected
    if isempty(parents)
        error('No parents selected for mutation.');
    end

    % Debugging: Display parent indices
    disp('Parent Indices:');
    disp(parents);

    % Extract parent chromosome
    parentChromosome = thisPopulation(parents, :);

    % Debugging: Display parent chromosome before mutation
    disp('Parent Chromosome Before Mutation:');
    disp(parentChromosome);

    % Ensure chromosome has the expected length
    if size(parentChromosome, 2) ~= nvars
        error(['Parent chromosome does not match the expected length of ', num2str(nvars), '.']);
    end

    % Extract relative dates and apply mutation
    relativeDates = parentChromosome(2:params.maxFlybys + 1);
    mutationStep = randn(size(relativeDates)) * 100; % Apply small random changes
    mutatedRelativeDates = relativeDates + mutationStep;

    % Ensure relative dates remain positive
    mutatedRelativeDates(mutatedRelativeDates < 0) = rand * 100;

    % Adjust relative dates to stay within mission duration
    if sum(mutatedRelativeDates) >= synodicPeriod
        mutatedRelativeDates = mutatedRelativeDates / sum(mutatedRelativeDates) * (synodicPeriod - rand * 100);
    end

    % Update chromosome with mutated relative dates
    mutatedChromosome = parentChromosome;
    mutatedChromosome(2:params.maxFlybys + 1) = mutatedRelativeDates;

    % Debugging: Display mutated chromosome
    disp('Mutated Chromosome:');
    disp(mutatedChromosome);

    % Ensure flyby planets are integers
    flybyPlanets = parentChromosome(params.maxFlybys + 2:end);
    mutatedChromosome(params.maxFlybys + 2:end) = round(flybyPlanets);
end
