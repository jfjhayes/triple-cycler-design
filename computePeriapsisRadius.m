function periapsis_radius = computePeriapsisRadius(v_inf_in, v_inf_out, mu)
% -------------------------------------------------------------------------
% Function: computePeriapsisRadius
% Purpose:
% Computes the periapsis radius of a gravity assist trajectory based on the 
% incoming and outgoing hyperbolic excess velocities and the planet's 
% gravitational parameter.
%
% Inputs:
%   v_inf_in  - Incoming hyperbolic excess velocity vector (km/s).
%   v_inf_out - Outgoing hyperbolic excess velocity vector (km/s).
%   mu        - Gravitational parameter of the flyby planet (km³/s²).
%
% Outputs:
%   periapsis_radius - Periapsis radius of the hyperbolic trajectory (km).
%
% Dependencies:
%   None.
% -------------------------------------------------------------------------

    % Compute the turning angle (delta) between the incoming and outgoing velocities
    delta = acos(dot(v_inf_in, v_inf_out) / (norm(v_inf_in) * norm(v_inf_out)));

    % Compute the magnitude of the incoming hyperbolic excess velocity
    v_inf_mag = norm(v_inf_in);

    % Compute the semi-major axis of the hyperbolic trajectory
    semi_major_axis = -mu / v_inf_mag^2;

    % Compute the periapsis radius of the hyperbolic trajectory
    periapsis_radius = semi_major_axis * (1 - 1 / cos(delta / 2));
end
