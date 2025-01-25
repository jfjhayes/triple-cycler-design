function planet = planet_properties(planet_id)
% -------------------------------------------------------------------------
% Function: planet_properties
% Purpose:
% This function returns the physical and orbital properties of a planet 
% needed for gravity assist computations. The planets are indexed as:
% Mercury (1) through Pluto (9).
%
% Inputs:
%   planet_id - Integer representing the planet ID (1–9)
%               Mercury (1), Venus (2), Earth (3), Mars (4), Jupiter (5),
%               Saturn (6), Uranus (7), Neptune (8), Pluto (9).
%
% Outputs:
%   planet - Structure containing the following properties:
%       - id: Planet ID
%       - mu: Gravitational parameter (km³/s²)
%       - radius: Mean radius (km)
%       - mass: Mass (kg)
%       - a: Semi-major axis of orbit (km)
%       - soi_radius: Sphere of influence radius (km)
%
% Dependencies:
%   compute_soi_radius - Subfunction to calculate the sphere of influence.
% -------------------------------------------------------------------------

    % Gravitational parameters and mass of the Sun
    sun.mu = 1.32712440018e11;  % Sun's gravitational parameter (km³/s²)
    sun.mass = 1.989e30;        % Sun's mass (kg)

    % Assign planet properties based on planet_id
    planet.id = planet_id;
    switch planet_id
        case 1 % Mercury
            planet.mu = 2.2032e4;  % km³/s²
            planet.radius = 2440;  % km
            planet.mass = 3.302e23; % kg
            planet.a = 0.387098 * 149597870.7;  % Semi-major axis (AU to km)
        case 2 % Venus
            planet.mu = 3.257e5;
            planet.radius = 6052;
            planet.mass = 4.8685e24;
            planet.a = 0.723332 * 149597870.7;
        case 3 % Earth
            planet.mu = 3.986004418e5;
            planet.radius = 6378;
            planet.mass = 5.97237e24;
            planet.a = 1.000000 * 149597870.7;
        case 4 % Mars
            planet.mu = 4.282837e4;
            planet.radius = 3396;
            planet.mass = 6.4171e23;
            planet.a = 1.523679 * 149597870.7;
        otherwise
            error('Planet ID not supported. Please provide an ID between 1 and 5.');
    end

    % Calculate the sphere of influence (SOI) radius
    planet.soi_radius = compute_soi_radius(planet, sun);

end

function r_soi = compute_soi_radius(planet, sun)
% -------------------------------------------------------------------------
% Subfunction: compute_soi_radius
% Purpose:
% Computes the sphere of influence (SOI) radius for a planet using the 
% formula: r_soi = a * (m_planet / m_sun)^(2/5).
%
% Inputs:
%   planet - Structure containing the planet's semi-major axis and mass.
%   sun    - Structure containing the Sun's mass.
%
% Outputs:
%   r_soi - Sphere of influence radius (km).
% -------------------------------------------------------------------------

    r_soi = planet.a * (planet.mass / sun.mass)^(2 / 5);
end
