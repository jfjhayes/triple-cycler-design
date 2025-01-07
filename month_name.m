function name = month_name(month_id)
    % Define the month names based on the month ID
    month_list = {'January', 'February', 'March', 'April', 'May', 'June', ...
                  'July', 'August', 'September', 'October', 'November', 'December'};
    
    % Ensure the month_id is within the valid range
    if month_id >= 1 && month_id <= length(month_list)
        name = month_list{month_id}; % Get the month name
    else
        error('Invalid month ID');
    end
end
