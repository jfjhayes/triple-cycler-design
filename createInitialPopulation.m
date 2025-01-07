function pop = createInitialPopulation(GenomeLength, ~, options)
    pop_size = options.PopulationSize;
    pop = zeros(pop_size, GenomeLength);
    num_dates = 6;

    % Constants
    date_min = J0(2010, 1, 1);
    date_max = J0(2050, 1, 1);
    min_days = 100;
    max_days = 500;

    for i = 1:pop_size
        valid = false;
        attempts = 0;
        max_attempts = 2000; % Reduced for faster iteration during testing

        while ~valid && attempts < max_attempts
            % Create first date
            dates = zeros(1, num_dates);
            dates(1) = date_min + rand() * (date_max - date_min - num_dates * min_days);

            % Add subsequent dates with spacing
            for j = 2:num_dates
                gap = min_days + rand() * (max_days - min_days);
                dates(j) = dates(j - 1) + gap;
            end

            % B-plane parameters
            b_mag = 1e3 + rand(1, 5) * (1e5 - 1e3);
            b_ang = -pi + rand(1, 5) * 2 * pi;

            % Construct chromosome
            chromosome = [dates reshape([b_mag; b_ang], 1, [])];

            % Debugging output
            fprintf('Generated dates: ');
            fprintf('%.2f ', dates);
            fprintf('\nDiffs: ');
            fprintf('%.2f ', diff(dates));
            fprintf('\n');

            % Feasibility checks directly embedded in generation logic
            if all(dates >= date_min) && all(dates <= date_max) && ...
               all(diff(dates) >= min_days) && all(diff(dates) <= max_days)
                pop(i, :) = chromosome;
                valid = true;
            end

            attempts = attempts + 1;
        end

        if ~valid
            warning('Failed to generate feasible individual after %d attempts', max_attempts);
            % Assign the last attempted chromosome to avoid gaps
            pop(i, :) = chromosome;
        end
    end
end

% Removed `isBasicallyFeasible` function since its checks are now inline


% function feasible = isBasicallyFeasible(chromosome, date_min, date_max, min_days, max_days)
%     dates = chromosome(1:6);
% 
%     feasible = true;
%     if any(dates < date_min) || any(dates > date_max)
%         feasible = false;
%         return;
%     end
% 
%     diffs = diff(dates);
%     if any(diffs <= 0) || any(diffs < min_days/365.25) || any(diffs > max_days/365.25)
%         feasible = false;
%         return;
%     end
% end