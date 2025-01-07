function [cost] = fitnessFunction(x, constants, constraints)
    try
        % Extract chromosome parameters
        [dates, b_params] = decodeChromosome(x);
        
        % Initialize costs and metrics
        taxi_dv = 0;
        penalty = 0;
        total_tof = 0;
        max_vinf = 0;
        
        % Process each trajectory segment
        for i = 1:size(dates,1)-1
            seg_params = struct(...
                'planet1', getPlanetID(i), ...
                'planet2', getPlanetID(i+1), ...
                'mu_sun', constants.mu_sun);
            
            [dv, tof, leg_data] = propagateSegment(...
                dates(i:i+1), b_params(i,1), b_params(i,2), seg_params);
            
            % Accumulate metrics
            taxi_dv = taxi_dv + dv;
            total_tof = total_tof + tof;
            vinf = norm(leg_data.v_in - leg_data.v_planet);
            max_vinf = max(max_vinf, vinf);
            
            % Calculate penalties
            if leg_data.altitude < constraints.alt_min
                penalty = penalty + 1e6*(constraints.alt_min - leg_data.altitude)^2;
            end
            if vinf > constraints.vinf_max
                penalty = penalty + 1e6*(vinf - constraints.vinf_max)^2;
            end
        end
        
        % Period closure constraint
        target_period = constants.synodic_period;
        actual_period = (dates(end) - dates(1))*86400;
        period_error = abs(actual_period - target_period);
        penalty = penalty + 1e6*(period_error/target_period)^2;
        
        % Final cost
        cost = taxi_dv + penalty;
        
        % Detailed output
        fprintf('\nTrajectory Evaluation:\n');
        fprintf('Total ΔV: %.2f km/s\n', taxi_dv);
        fprintf('Max v∞: %.2f km/s\n', max_vinf);
        fprintf('Total ToF: %.1f days\n', total_tof/86400);
        fprintf('Period error: %.1f%%\n', 100*period_error/target_period);
        fprintf('Penalty: %.2e\n', penalty);
        fprintf('Final cost: %.2e\n\n', cost);
        
    catch ME
        fprintf('Fitness evaluation failed: %s\n', ME.message);
        cost = 1e10;
    end
end

function penalty = evaluateConstraints(leg_data, constraints)
    penalty = 0;
    
    % Altitude constraints
    if leg_data.altitude < constraints.alt_min
        penalty = penalty + 1e6*(constraints.alt_min - leg_data.altitude)^2;
    elseif leg_data.altitude > constraints.alt_max
        penalty = penalty + 1e6*(leg_data.altitude - constraints.alt_max)^2;
    end
    
    % V-infinity constraints
    vinf = norm(leg_data.v_in - leg_data.v_planet);
    if vinf < constraints.vinf_min
        penalty = penalty + 1e6*(constraints.vinf_min - vinf)^2;
    elseif vinf > constraints.vinf_max
        penalty = penalty + 1e6*(vinf - constraints.vinf_max)^2;
    end
    
    % Transfer time constraints
    if leg_data.tof < constraints.tof_min
        penalty = penalty + 1e6*(constraints.tof_min - leg_data.tof)^2;
    elseif leg_data.tof > constraints.tof_max
        penalty = penalty + 1e6*(leg_data.tof - constraints.tof_max)^2;
    end
end

function error = evaluateClosureError(date1, date2, constants)
    % Check if trajectory closes near synodic period
    actual_period = (datenum(date2) - datenum(date1))*86400;
    error = abs(actual_period - constants.synodic_period)/constants.synodic_period;
end