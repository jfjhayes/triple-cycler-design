%% Genetic Algorithm Implementation for Optimising Trajectories
% -------------------------------------------------------------------------
% Purpose:
% This script implements a Genetic Algorithm (GA) to optimise planetary
% trajectories for a low-excess-speed triple cycler trajectory. The GA is
% configured with custom constraints, mutation, and output functions.
%
% Workflow:
% 1. Define mission parameters and generate an initial population.
% 2. Configure GA options and fitness function.
% 3. Run the GA to find the optimal trajectory chromosome.
% 4. Decode the best chromosome and compute the trajectory for analysis.
%
% Inputs:
% - Mission parameters defined in `defineParams`.
%
% Outputs:
% - Best trajectory chromosome and associated fitness value (Delta-V).
% - Detailed trajectory logs for further analysis.
% -------------------------------------------------------------------------

% Define Parameters
params = defineParams();

% Generate Initial Population
popSize = 1000; % Population size
population = createInitialPopulation(popSize, params);

% Define Chromosome Structure
numVars = params.maxFlybys * 2 + 1; % Start date + flyby dates + flyby planets

% Define Fitness Function
fitnessFcn = @(chromosome) fitnessFunction(chromosome, params);

% Configure GA Options
options = optimoptions('ga', ...
    'MutationFcn', @(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation) ...
                  customMutation(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation, params), ...
    'PopulationSize', popSize, ...
    'MaxGenerations', 500, ...
    'MaxStallGenerations', 50, ...
    'FunctionTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-3, ...
    'CrossoverFraction', 0.9, ...
    'OutputFcn', @(options, state, flag) gaOutputFunction(options, state, flag, params), ...
    'PlotFcn', {@gaplotbestf, @gaplotstopping, @gaplotscorediversity}, ...
    'Display', 'iter', ...
    'UseParallel', true);

% Define Bounds
lb = [params.start_date, zeros(1, params.maxFlybys), 2 * ones(1, params.maxFlybys)]; % Lower bounds
ub = [params.end_date, ones(1, params.maxFlybys) * (params.end_date - params.start_date), 4 * ones(1, params.maxFlybys)]; % Upper bounds

% Define Integer Constraints
intcon = params.maxFlybys + 2 : numVars; % Flyby planets are integers

% Preallocate Storage for Logging
saveInterval = 10; % Save logs every 10 generations
populationLog = struct('generation', [], 'chromosomes', [], 'fitness', []);

% Run Genetic Algorithm
disp('Running Genetic Algorithm...');
[bestChromosome, bestFitness, ~, output] = ga(fitnessFcn, numVars, [], [], [], [], lb, ub, [], intcon, options);

% Load the saved population log (for debugging or analysis)
load('final_population_log_Mars.mat', 'populationLog');

% Extract the best chromosome and fitness from the last generation
finalGen = length(populationLog);
[bestFitness, bestIndex] = min(populationLog(finalGen).fitness);
bestChromosome = populationLog(finalGen).chromosomes(bestIndex, :);

% Decode the best chromosome
bestStartDate = bestChromosome(1);
bestRelativeDates = bestChromosome(2:params.maxFlybys+1);
bestFlybyPlanets = bestChromosome(params.maxFlybys+2:end);

% Compute trajectory for the best solution
[legs, totalDeltaV, isFeasible, logData] = computeTrajectory(bestStartDate, bestRelativeDates, bestFlybyPlanets, params);

% Display Results
disp('Best Chromosome:');
disp(bestChromosome);
disp(['Best Fitness (Delta-V): ', num2str(bestFitness)]);
disp(['Start Date: ', datestr(datetime(bestStartDate, 'ConvertFrom', 'juliandate'), 'yyyy-mm-dd HH:MM:SS')]);
disp(['Relative Flyby Dates: ', num2str(bestRelativeDates)]);
disp(['Flyby Planets: ', num2str(bestFlybyPlanets)]);
disp('Detailed Trajectory Log:');
disp(logData);
