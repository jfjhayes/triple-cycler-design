function [dates, b_params] = decodeChromosome(x)
    % Constants  
    num_dates = 6;
    min_transfer_days = 100;
    max_transfer_days = 500;
    date_min = J0(2010, 1, 1);
    date_max = J0(2050, 1, 1);
    
    % Extract dates
    dates = x(1:num_dates)';
    
    % Extract B-plane parameters
    b_params = reshape(x(num_dates+1:num_dates+10), 2, 5)';
    
    % Return actual dates but let constraints handle validation
    % This allows GA to evolve invalid solutions with penalties
end