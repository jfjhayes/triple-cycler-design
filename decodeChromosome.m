function [start_date, flyby_dates, flyby_planets, end_date] = decodeChromosome(chromosome, params)
% -------------------------------------------------------------------------
% Function: decodeChromosome
% Purpose:
% This function decodes a chromosome representing a trajectory into its 
% components: start date, relative flyby dates, flyby planets, and end date.
% The decoding respects mission constraints provided in the `params` structure.
%
% Inputs:
%   chromosome - A vector encoding the trajectory parameters:
%                [start_date, relative_flyby_dates, flyby_planets].
%   params     - Structure containing mission-specific parameters:
%                - maxFlybys: Maximum number of flybys allowed in the mission.
%                - synodicPeriod: Duration of the mission in days.
%
% Outputs:
%   start_date    - The Julian start date of the mission.
%   flyby_dates   - Relative dates of the flybys (days from `start_date`).
%   flyby_planets - Sequence of planets encountered during flybys.
%   end_date      - The fixed end date of the mission (start_date + synodicPeriod).
% -------------------------------------------------------------------------

    % Extract the start date (Julian)
    start_date = chromosome(1);

    % Extract relative flyby dates (days from the start date)
    flyby_dates = chromosome(2:params.maxFlybys+1);

    % Extract the sequence of flyby planets
    flyby_planets = chromosome(params.maxFlybys+2:end);

    % Calculate the fixed end date based on the synodic period
    end_date = start_date + params.synodicPeriod;
end
