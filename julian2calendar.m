function [year, month, day] = julian2calendar(jd)
    % Converts Julian Date to calendar date using inverse of J0 algorithm
    % Valid for years 1901-2099
    
    jd = jd - 1721013.5;  % Remove J0 offset
    
    % First approximation
    year = floor(jd/365.25);
    
    % Refine year
    while J0(year+1, 1, 1) - 1721013.5 <= jd
        year = year + 1;
    end
    
    % Find month
    month = 1;
    while month < 12 && J0(year, month+1, 1) - 1721013.5 <= jd
        month = month + 1;
    end
    
    % Calculate day
    day = floor(jd - (J0(year, month, 1) - 1721013.5) + 1);

    % Validate results
    validateattributes(year, {'numeric'}, {'integer', '>=', 1901, '<=', 2099});
    validateattributes(month, {'numeric'}, {'integer', '>=', 1, '<=', 12});
    validateattributes(day, {'numeric'}, {'integer', '>=', 1, '<=', 31});
end