function E = kepler_E(e, M)
% -------------------------------------------------------------------------
% Function: kepler_E
% Purpose:
% Solves Kepler's equation:
%   E - e*sin(E) = M
% for the eccentric anomaly (E) using Newton's method. The solution 
% iterates until the desired error tolerance is achieved.
%
% Inputs:
%   e - Eccentricity (unitless).
%   M - Mean anomaly (radians).
%
% Outputs:
%   E - Eccentric anomaly (radians).
%
% Dependencies:
%   None.
% -------------------------------------------------------------------------

    % Set the error tolerance for convergence
    tolerance = 1.e-8;

    % Initial guess for E based on the mean anomaly (M)
    if M < pi
        E = M + e / 2;  % For M < π, add e/2
    else
        E = M - e / 2;  % For M ≥ π, subtract e/2
    end

    % Iterate using Newton's method until the solution converges
    ratio = 1;  % Initialise ratio for the loop
    while abs(ratio) > tolerance
        % Compute the difference and update E
        ratio = (E - e * sin(E) - M) / (1 - e * cos(E));
        E = E - ratio;
    end
end
