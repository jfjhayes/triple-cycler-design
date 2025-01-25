function coe = coe_from_sv(R, V, mu)
% -------------------------------------------------------------------------
% Function: coe_from_sv
% Purpose:
% Computes the classical orbital elements (COEs) from a given state vector 
% (position and velocity) using Algorithm 4.1.
%
% Inputs:
%   R  - Position vector in the geocentric equatorial frame (km).
%   V  - Velocity vector in the geocentric equatorial frame (km/s).
%   mu - Gravitational parameter of the central body (km³/s²).
%
% Outputs:
%   coe - Vector of classical orbital elements:
%         [h, e, RA, incl, w, TA, a], where:
%         - h: Specific angular momentum (km²/s)
%         - e: Eccentricity
%         - RA: Right Ascension of the Ascending Node (radians)
%         - incl: Inclination (radians)
%         - w: Argument of Periapsis (radians)
%         - TA: True Anomaly (radians)
%         - a: Semi-major Axis (km)
%
% Dependencies:
%   None.
% -------------------------------------------------------------------------

    % Define a small numerical tolerance for eccentricity checks
    eps = 1.e-10;

    % Compute magnitudes of position and velocity vectors
    r = norm(R);  % Magnitude of the position vector (km)
    v = norm(V);  % Magnitude of the velocity vector (km/s)

    % Compute the radial velocity component
    vr = dot(R, V) / r;

    % Compute the angular momentum vector and its magnitude
    H = cross(R, V);  % Angular momentum vector (km²/s)
    h = norm(H);      % Specific angular momentum (km²/s)

    % Compute the inclination (Equation 4.7)
    incl = acos(H(3) / h);

    % Compute the node line vector and its magnitude (Equation 4.8)
    N = cross([0, 0, 1], H);  % Node line vector (km²/s)
    n = norm(N);

    % Compute the Right Ascension of the Ascending Node (RA) (Equation 4.9)
    if n ~= 0
        RA = acos(N(1) / n);
        if N(2) < 0
            RA = 2 * pi - RA;  % Adjust for the quadrant
        end
    else
        RA = 0;  % Undefined if node line magnitude is zero
    end

    % Compute the eccentricity vector and its magnitude (Equation 4.10)
    E = (1 / mu) * ((v^2 - mu / r) * R - r * vr * V);
    e = norm(E);

    % Compute the Argument of Periapsis (w) (Equation 4.12)
    if n ~= 0
        if e > eps
            w = acos(dot(N, E) / (n * e));
            if E(3) < 0
                w = 2 * pi - w;  % Adjust for the quadrant
            end
        else
            w = 0;  % If eccentricity is zero, w is undefined
        end
    else
        w = 0;  % If node line magnitude is zero, w is undefined
    end

    % Compute the True Anomaly (TA) (Equation 4.13a)
    if e > eps
        TA = acos(dot(E, R) / (e * r));
        if vr < 0
            TA = 2 * pi - TA;  % Adjust for the quadrant
        end
    else
        cp = cross(N, R);  % Cross product to determine direction
        if cp(3) >= 0
            TA = acos(dot(N, R) / (n * r));
        else
            TA = 2 * pi - acos(dot(N, R) / (n * r));
        end
    end

    % Compute the Semi-major Axis (a) (Equation 4.62)
    % a < 0 for hyperbolic orbits
    a = h^2 / mu / (1 - e^2);

    % Store the computed orbital elements in the output vector
    coe = [h, e, RA, incl, w, TA, a];
end
