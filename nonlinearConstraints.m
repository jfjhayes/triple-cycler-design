function [c, ceq] = nonlinearConstraints(x, constants, constraints)
   % Initialize constraints
   num_segments = 5;  % Number of trajectory segments
   c = zeros(4*num_segments, 1);    % 4 inequality constraints per segment
   ceq = zeros(2*(num_segments-1) + 1, 1);  % 2 continuity per junction + period

   try
       % Decode chromosome
       [dates, b_params] = decodeChromosome(x);

       % In nonlinearConstraints:
       fprintf('Decoded dates:\n');
       disp(dates);
       fprintf('Date differences:\n');
       disp(diff(dates));
       
       % Time constraints 
       diffs = diff(dates)*86400;
        c_time = [
            constraints.tof_min - diffs;  % Minimum time
            diffs - constraints.tof_max   % Maximum time
        ];
       
       c_idx = 1;
       ceq_idx = 1;
       prev_leg = [];
       
       for i = 1:size(dates,1)-1
           seg_params = struct(...
               'planet1', getPlanetID(i), ...
               'planet2', getPlanetID(i+1), ...
               'mu_sun', constants.mu_sun);
           
           [~, ~, leg_data] = propagateSegment(...
               dates(i:i+1), b_params(i,1), b_params(i,2), seg_params);
           
           if ~isempty(prev_leg)
               % Continuity constraints
               ceq(ceq_idx:ceq_idx+1) = [...
                   norm(leg_data.r_planet - prev_leg.r_target) - constraints.pos_tol;
                   norm(leg_data.v_planet - prev_leg.v_target) - constraints.vel_tol];
               ceq_idx = ceq_idx + 2;
           end
           
           % Inequality constraints
           c(c_idx:c_idx+3) = [...
               leg_data.altitude - constraints.alt_max;
               constraints.alt_min - leg_data.altitude;
               norm(leg_data.v_in - leg_data.v_planet) - constraints.vinf_max;
               norm(leg_data.v_out - leg_data.v_planet) - constraints.vinf_max];
           
           c_idx = c_idx + 4;
           prev_leg = leg_data;
       end
       
       % Period constraint
       period_error = (dates(end) - dates(1))*86400 - constants.synodic_period;
       ceq(end) = period_error;
       
       % Append time constraints
       c = [c; c_time];
       
   catch ME
       fprintf('Error in constraints: %s\n', ME.message);
       fprintf('Error stack:\n');
       for k = 1:length(ME.stack)
           fprintf('  %s: line %d\n', ME.stack(k).name, ME.stack(k).line);
       end
       c(:) = Inf;
       ceq(:) = Inf;
   end
   
   % Validate outputs
   if any(isnan(c)) || any(isnan(ceq)) || ~isreal(c) || ~isreal(ceq)
       c(:) = Inf;
       ceq(:) = Inf;
   end
end