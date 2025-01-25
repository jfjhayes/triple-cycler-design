function [r, v] = sv_from_coe(coe, mu)
% -------------------------------------------------------------------------
% Function: sv_from_coe
% Purpose:
% Computes the state vector (position and velocity) in the geocentric 
% equatorial frame from the classical orbital elements (COEs).
%
% Inputs:
%   coe - Vector of classical orbital elements:
%         [h, e, RA, incl, w, TA], where:
%         - h: Specific angular momentum (km²/s)
%         - e: Eccentricity
%         - RA: Right Ascension of the Ascending Node (radians)
%         - incl: Inclination (radians)
%         - w: Argument of Perigee (radians)
%         - TA: True Anomaly (radians)
%   mu  - Gravitational parameter of the central body (km³/s²).
%
% Outputs:
%   r - Position vector in the geocentric equatorial frame (km).
%   v - Velocity vector in the geocentric equatorial frame (km/s).
%
% Dependencies:
%   None.
% -------------------------------------------------------------------------

    % Extract classical orbital elements
    h = coe(1);    % Specific angular momentum (km²/s)
    e = coe(2);    % Eccentricity
    RA = coe(3);   % Right Ascension of the Ascending Node (radians)
    incl = coe(4); % Inclination (radians)
    w = coe(5);    % Argument of Perigee (radians)
    TA = coe(6);   % True Anomaly (radians)

    % Compute position and velocity in the perifocal frame (Equations 4.45 and 4.46)
    rp = (h^2 / mu) * (1 / (1 + e * cos(TA))) * ...
         (cos(TA) * [1; 0; 0] + sin(TA) * [0; 1; 0]);
    vp = (mu / h) * (-sin(TA) * [1; 0; 0] + (e + cos(TA)) * [0; 1; 0]);

    % Rotation matrix to transform from perifocal to geocentric equatorial frame
    % Rotation about the z-axis through the angle RA
    R3_W = [ cos(RA),  sin(RA), 0;
            -sin(RA),  cos(RA), 0;
                 0,         0, 1];

    % Rotation about the x-axis through the angle incl
    R1_i = [1,        0,         0;
            0,  cos(incl),  sin(incl);
            0, -sin(incl),  cos(incl)];

    % Rotation about the z-axis through the angle w
    R3_w = [ cos(w),  sin(w), 0;
            -sin(w),  cos(w), 0;
                 0,       0, 1];

    % Combined transformation matrix (Equation 4.49)
    Q_pX = (R3_w * R1_i * R3_W)';

    % Transform position and velocity from perifocal to geocentric equatorial frame
    r = Q_pX * rp;  % Position vector (column vector)
    v = Q_pX * vp;  % Velocity vector (column vector)

    % Convert position and velocity to row vectors
    r = r';
    v = v';
end
