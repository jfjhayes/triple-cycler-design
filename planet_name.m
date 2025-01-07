function name = planet_name(planet_id)
    % Define the planet names based on the ID
    planet_list = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'};
    
    % Ensure the planet_id is within the valid range
    if planet_id >= 1 && planet_id <= length(planet_list)
        name = planet_list{planet_id}; % Get the planet name
    else
        error('Invalid planet ID');
    end
end