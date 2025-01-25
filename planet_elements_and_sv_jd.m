function [coe, r, v, jd] = planet_elements_and_sv_jd(planet_id, jd)
% -------------------------------------------------------------------------
% Function: planet_elements_and_sv_jd
% Purpose:
% Computes the classical orbital elements and the heliocentric state vector 
% (position and velocity) of a planet for a given Julian date.
%
% Inputs:
%   planet_id - Integer ID representing the planet (1–9 for Mercury–Pluto).
%               For example: 3 = Earth, 4 = Mars.
%   jd        - Julian date for which the computation is performed.
%
% Outputs:
%   coe - Orbital elements of the planet:
%         [h, e, RA, incl, w, TA, a, w_hat, L, M, E]
%         where:
%         - h: Specific angular momentum (km²/s)
%         - e: Eccentricity
%         - RA: Right Ascension of the Ascending Node (degrees)
%         - incl: Inclination (degrees)
%         - w: Argument of Periapsis (degrees)
%         - TA: True Anomaly (degrees)
%         - a: Semi-major Axis (km)
%         - w_hat: Longitude of Perihelion (degrees)
%         - L: Mean Longitude (degrees)
%         - M: Mean Anomaly (degrees)
%         - E: Eccentric Anomaly (degrees)
%   r   - Heliocentric position vector in Cartesian coordinates (km).
%   v   - Heliocentric velocity vector in Cartesian coordinates (km/s).
%   jd  - Julian date (for consistency).
%
% Dependencies:
%   - planetary_elements (subfunction): Provides J2000 orbital elements and rates.
%   - zero_to_360 (subfunction): Reduces angular values to the 0–360 degree range.
%   - sv_from_coe: Converts orbital elements to state vectors.
%   - kepler_E: Solves Kepler's equation for eccentric anomaly.
% -------------------------------------------------------------------------

    % Gravitational parameter of the Sun (km³/s²)
    mu = 1.32712440018e11;

    % Conversion factor: Degrees to Radians
    deg = pi / 180;

    % Julian centuries since J2000 epoch
    t0 = (jd - 2451545) / 36525;

    % Get J2000 orbital elements and their rates of change
    [J2000_coe, rates] = planetary_elements(planet_id);

    % Compute the orbital elements at the given Julian date
    elements = J2000_coe + rates * t0;

    % Extract semi-major axis and eccentricity
    a = elements(1);  % Semi-major axis (km)
    e = elements(2);  % Eccentricity

    % Compute angular momentum
    h = sqrt(mu * a * (1 - e^2));

    % Reduce angular elements to the 0–360 degree range
    incl = elements(3);  % Inclination
    RA = zero_to_360(elements(4));  % Right Ascension
    w_hat = zero_to_360(elements(5));  % Longitude of Perihelion
    L = zero_to_360(elements(6));  % Mean Longitude
    w = zero_to_360(w_hat - RA);  % Argument of Periapsis
    M = zero_to_360(L - w_hat);  % Mean Anomaly

    % Solve Kepler's Equation to find Eccentric Anomaly (E)
    E = kepler_E(e, M * deg);  % Eccentric anomaly (radians)

    % Compute True Anomaly (TA)
    TA = zero_to_360(2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) / deg);

    % Store orbital elements in the output
    coe = [h, e, RA, incl, w, TA, a, w_hat, L, M, E / deg];

    % Convert orbital elements to state vectors
    [r, v] = sv_from_coe([h, e, RA * deg, incl * deg, w * deg, TA * deg], mu);
end

function [J2000_coe, rates] = planetary_elements(planet_id)
% -------------------------------------------------------------------------
% Subfunction: planetary_elements
% Purpose:
% Retrieves the J2000 epoch orbital elements and their centennial rates 
% of change for the specified planet.
%
% Inputs:
%   planet_id - Integer representing the planet (1–9 for Mercury–Pluto).
%
% Outputs:
%   J2000_coe - Vector of J2000 epoch orbital elements:
%               [a, e, incl, RA, w_hat, L] (in km, degrees).
%   rates     - Vector of rates of change per Julian century:
%               [a_dot, e_dot, incl_dot, RA_dot, w_hat_dot, L_dot].
% -------------------------------------------------------------------------

    % Orbital elements and rates for the nine planets
    J2000_elements = ...
        [ 0.38709893 0.20563069 7.00487 48.33167 77.45645 252.25084;
          0.72333199 0.00677323 3.39471 76.68069 131.53298 181.97973;
          1.00000011 0.01671022 0.00005 -11.26064 102.94719 100.46435;
          1.52366231 0.09341233 1.85061 49.57854 336.04084 355.45332;
          5.20336301 0.04839266 1.30530 100.55615 14.75385 34.40438;
          9.53707032 0.05415060 2.48446 113.71504 92.43194 49.94432;
          19.19126393 0.04716771 0.76986 74.22988 170.96424 313.23218;
          30.06896348 0.00858587 1.76917 131.72169 44.97135 304.88003;
          39.48168677 0.24880766 17.14175 110.30347 224.06676 238.92881 ];

    cent_rates = ...
        [ 0.00000066 0.00002527 -23.51 -446.30 573.57 538101628.29;
          0.00000092 -0.00004938 -2.86 -996.89 -108.80 210664136.06;
         -0.00000005 -0.00003804 -46.94 -18228.25 1198.28 129597740.63;
         -0.00007221 0.00011902 -25.47 -1020.19 1560.78 68905103.78;
          0.00060737 -0.00012880 -4.15 1217.17 839.93 10925078.35;
         -0.00301530 -0.00036762 6.11 -1591.05 -1948.89 4401052.95;
          0.00152025 -0.00019150 -2.09 -1681.4 1312.56 1542547.79;
         -0.00125196 0.00002514 -3.64 -151.25 -844.43 786449.21;
         -0.00076912 0.00006465 11.07 -37.33 -132.25 522747.90 ];

    % Extract data for the specified planet
    J2000_coe = J2000_elements(planet_id, :);
    rates = cent_rates(planet_id, :);

    % Convert semi-major axis from AU to km
    au = 149597871;
    J2000_coe(1) = J2000_coe(1) * au;
    rates(1) = rates(1) * au;

    % Convert arcseconds to degrees
    rates(3:6) = rates(3:6) / 3600;
end

function y = zero_to_360(x)
% -------------------------------------------------------------------------
% Subfunction: zero_to_360
% Purpose:
% Reduces an angle to lie within the range 0–360 degrees.
%
% Inputs:
%   x - Input angle in degrees.
%
% Outputs:
%   y - Normalised angle in the range 0–360 degrees.
% -------------------------------------------------------------------------

    if x >= 360
        x = x - fix(x / 360) * 360;
    elseif x < 0
        x = x - (fix(x / 360) - 1) * 360;
    end
    y = x;
end
