function [v_out, alt, turn] = computeGravityAssist(v_in, v_planet, planet_data, b_mag, b_angle)
    % COMPUTEGRAVITYASSIST Computes gravity assist using B-plane targeting
    %
    % Uses formulation from Battin's "Astrodynamics" for accurate periapsis
    % and turn angle calculations in the B-plane reference frame.
    %
    % Inputs:
    %   v_in        - Incoming velocity vector [km/s]
    %   v_planet    - Planet velocity vector [km/s]  
    %   planet_data - Structure with fields:
    %                 .mu      - Gravitational parameter [km^3/s^2]
    %                 .radius  - Planet radius [km]
    %   b_mag       - B-plane vector magnitude [km]
    %   b_angle     - B-plane vector angle [rad]
    %
    % Outputs:
    %   v_out       - Post-flyby velocity vector [km/s]
    %   alt         - Flyby altitude [km]
    %   turn        - Turn angle achieved [rad]

    % Input validation 
    validateattributes(v_in, {'numeric'}, {'vector', 'numel', 3, 'finite', 'real'});
    validateattributes(v_planet, {'numeric'}, {'vector', 'numel', 3, 'finite', 'real'});
    validateattributes(b_mag, {'numeric'}, {'scalar', '>', 0, 'finite'});
    validateattributes(b_angle, {'numeric'}, {'scalar', 'finite', 'real'});

    % Ensure column vectors
    v_in = v_in(:);
    v_planet = v_planet(:);

    % Calculate v-infinity
    v_inf_in = v_in - v_planet;
    v_inf_mag = norm(v_inf_in);

    % Check for physical validity
    if v_inf_mag < sqrt(planet_data.mu/planet_data.radius)
        error('V-infinity below escape velocity - capture would occur');
    end
    
    if v_inf_mag < 1e-10
        error('Zero relative velocity between spacecraft and planet');
    end

    % Construct B-plane basis
    s_hat = v_inf_in/v_inf_mag;
    
    % Choose reference vector avoiding parallel alignment
    if abs(dot(s_hat, [0;0;1])) > 0.9
        k = [1;0;0];
    else
        k = [0;0;1];
    end
    
    t_hat = cross(s_hat, k);
    t_hat = t_hat/norm(t_hat);
    r_hat = cross(s_hat, t_hat);

    % Compute B-vector
    B = b_mag * (cos(b_angle)*t_hat + sin(b_angle)*r_hat);
    B_hat = B/norm(B);

    % Calculate periapsis and turn angle using Battin's formulation
    e = sqrt(1 + (b_mag*v_inf_mag^2/planet_data.mu)^2);  % Hyperbolic eccentricity
    rp = planet_data.mu/(v_inf_mag^2) * (e - 1);         % Periapsis radius
    turn = 2 * asin(1/e);                                 % Turn angle

    % Calculate outbound v-infinity using Rodrigues formula
    v_inf_out = v_inf_mag * (s_hat*cos(turn)) + ...
                cross(B_hat, s_hat)*sin(turn) + ...
                B_hat*dot(B_hat, s_hat)*(1-cos(turn));

    % Compute output parameters
    v_out = v_inf_out + v_planet;
    alt = rp - planet_data.radius;

    % Validate results
    if alt < 100  % Minimum safe altitude
        error('Computed altitude below minimum safe value');
    end
    
    if abs(norm(v_inf_out) - v_inf_mag) > 1e-10
        error('V-infinity magnitude not preserved');
    end
end