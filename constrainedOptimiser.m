%% Fmincon refining GA results
% -------------------------------------------------------------------------
% Purpose:
% This script refines a trajectory derived from the GA for a 
% low-excess-speed triple cycler between Venus, Earth, and Mars. The trajectory is optimised
% to minimise taxi Delta-V using `fmincon`, while adhering to time and planetary constraints.
%
% Inputs:
% - final_population_log_Earth.mat OR inal_population_log_Mars.mat: Contains the GA output (`populationLog`), including 
%   chromosome details for trajectory refinement.
% - Params from `defineParams.m`: Include mission-specific parameters like bounds and
%   synodic periods.
%
% Outputs:
% - refined_solution_Earth.mat: A MATLAB .mat file storing the refined trajectory details.
% - Console output: Displays initial and refined trajectory details for visual verification.
%
% Main Functions:
% - decodeChromosome: Decodes the GA chromosome into start date, flyby dates, and planet sequence.
% - computeObjective: Evaluates Delta-V for a given trajectory.
% - fmincon: Refines the trajectory's dates within defined bounds.
%
% Dependencies:
% - fmincon
% - Planetary state vectors and gravity assist calculations
% -------------------------------------------------------------------------

% Load parameters and GA output
params = defineParams();  % Mission-specific settings and constants
load('final_population_log_Mars.mat', 'populationLog');  % GA output

% Extract the final generation and the best trajectory chromosome
final_gen = populationLog(end);  % Access last generation in population log
[~, best_idx] = min(final_gen.fitness);  % Find the chromosome with the best fitness
best_chromosome = final_gen.chromosomes(best_idx, :);  % Extract best chromosome

% Decode chromosome to trajectory details
[start_date, flyby_dates, flyby_planets, end_date] = decodeChromosome(best_chromosome, params);

% Display the decoded initial trajectory for debugging
disp('Initial Trajectory:');
disp(['  Start Date (Julian): ', num2str(start_date)]);
disp(['  End Date (Julian): ', num2str(end_date)]);
disp('  Flyby Planets (GA guess):');
disp(flyby_planets);

% Define bounds for optimisation
perturbation_limit = 30;  % Limit of ±30 days for flyby dates
lb = max(0, flyby_dates - perturbation_limit);  % Ensure lower bounds stay positive
ub = min(params.synodicPeriod, flyby_dates + perturbation_limit);  % Upper bound within allowable mission duration

% Objective function to minimise Delta-V, including physical Delta-V and penalties
objective = @(flyby_dates) computeObjective(start_date, flyby_dates, flyby_planets, params);

% Configure options for the optimisation solver
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
                       'MaxFunctionEvaluations', 1e5, 'MaxIterations', 1e4, ...
                       'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-10);

% Optimise trajectory by refining flyby dates while keeping planetary sequence fixed
fprintf('Starting trajectory refinement...\n');
refined_flyby_dates = fmincon(objective, flyby_dates, [], [], [], [], lb, ub, [], options);

% Package the refined trajectory solution
refined_solution.start_date = start_date;
refined_solution.flyby_dates = refined_flyby_dates;
refined_solution.flyby_planets = flyby_planets;  % Retain GA-derived sequence
refined_solution.end_date = end_date;  % End date remains fixed

% Display and save the refined trajectory solution
disp('Refined Trajectory:');
disp(refined_solution);
save('refined_solution_Earth.mat', 'refined_solution');  % Save results
