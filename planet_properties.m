function planet = planet_properties(planet_id)
    % Returns physical properties for planets needed for gravity assists
    % Input: planet_id (1-9 corresponding to Mercury through Pluto)
    % Output: structure containing planetary properties
    
    % Define sun properties (needed for SOI calculation)
    sun.mu = 1.32712440018e11;
    sun.mass = 1.989e30;
    
    switch planet_id
        case 1 % Mercury
            planet.mu = 2.2032e4;
            planet.radius = 2440;
            planet.mass = 3.302e23;
            planet.a = 0.387098 * 149597870.7;  % AU to km
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
        case 5 % Jupiter
            planet.mu = 1.26686534e8;
            planet.radius = 71492;
            planet.mass = 1.8982e27;
            planet.a = 5.204267 * 149597870.7;
        otherwise
            error('Planet ID not supported');
    end
    
    % Add SOI calculation
    planet.soi_radius = compute_soi_radius(planet, sun);
end

function r_soi = compute_soi_radius(planet, sun)
    % Compute sphere of influence radius
    r_soi = planet.a * (planet.mass/sun.mass)^(2/5);
end