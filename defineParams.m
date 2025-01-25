function params = defineParams()
% -------------------------------------------------------------------------
% Function: defineParams
% Purpose:
% This function initialises a parameter structure containing general settings 
% for the Genetic Algorithm (GA), constrained optimiser, and trajectory computations. It includes time 
% bounds, flyby constraints, planetary data, and gravitational constants.
%
% Inputs:
% None.
%
% Outputs:
% params - Struct containing the following fields:
%   - start_date: Earliest allowable mission start date (Julian Date).
%   - end_date: Latest allowable mission end date (Julian Date).
%   - maxFlybys: Maximum number of flybys allowed in a chromosome.
%   - min_altitude: Minimum altitude for gravity assists (km).
%   - synodicPeriod: Twice the synodic period (days) for Venus-Earth-Mars alignments.
%   - planets: Array of structures containing planetary data (Venus, Earth, Mars).
%   - mu_sun: Gravitational parameter of the Sun (km^3/s^2).
%   - startPlanet: The initial planet for the mission (Earth by default).
% -------------------------------------------------------------------------


    % General mission parameters
    params.start_date = J0(2025, 1, 1); % Convert 1 January 2025 to Julian date
    params.end_date = J0(2032, 1, 1);   % Convert 1 January 2032 to Julian date
    params.maxFlybys = 4;               % Maximum flybys per trajectory
    params.min_altitude = 100;          % Minimum allowable flyby altitude (km)
    params.synodicPeriod = 365.25 * 6.4 * 2;  % Twice Venus-Earth-Mars synodic period

    % Load planetary properties for Venus, Earth, and Mars
    params.planets = [planet_properties(2), ... % Venus
                      planet_properties(3), ... % Earth
                      planet_properties(4)];    % Mars

    % Gravitational parameter of the Sun
    params.mu_sun = 1.32712440018e11;  % Sun's gravitational parameter (km^3/s^2)

    % Set the starting planet (Earth by default)
    params.startPlanet = 3; % Earth = 3, Mars = 4 (can be changed for inbound/outbound)
end
