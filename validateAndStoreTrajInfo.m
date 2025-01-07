function validateAndStoreTrajInfo(leg_data, v_dep, v_arr, r1, r2, v1, v2, mu_sun)
    % Compute orbital elements and v-infinity
    r1_norm = norm(r1);
    v_inf_dep = norm(v_dep - v1);
    v_inf_arr = norm(v_arr - v2);
    transfer_angle = acos(dot(r1,r2)/(r1_norm*norm(r2)));
    specific_energy = norm(v_dep)^2/2 - mu_sun/r1_norm;
    sma = -mu_sun/(2*specific_energy);
    
    % Store validation results
    if transfer_angle < pi
        transfer_type = 'Type I';
    else
        transfer_type = 'Type II';
    end
     
    leg_data.validation = struct(...
        'sma', sma, ...
        'specific_energy', specific_energy, ...
        'transfer_angle', transfer_angle, ...
        'transfer_type', transfer_type, ...
        'v_inf_dep', v_inf_dep, ...
        'v_inf_arr', v_inf_arr, ...
        'is_valid', v_inf_dep < 20 && v_inf_arr < 20 && specific_energy < 0 && sma > 0);
        
    % Output diagnostics  
    fprintf('\nTrajectory validation:\n');
    fprintf('SMA: %.2e km\n', sma);
    fprintf('V-infinity departure: %.2f km/s\n', v_inf_dep);
    fprintf('V-infinity arrival: %.2f km/s\n', v_inf_arr); 
    fprintf('Transfer angle: %.1f degrees\n', rad2deg(transfer_angle));
    fprintf('Transfer type: %s\n', leg_data.validation.transfer_type);
    
    if ~leg_data.validation.is_valid
        warning('Invalid trajectory detected');
        fprintf('Validation failures:\n');
        if v_inf_dep > 20 || v_inf_arr > 20
            fprintf('- Excessive v-infinity\n');
        end
        if specific_energy >= 0
            fprintf('- Hyperbolic escape trajectory\n'); 
        end
        if sma <= 0
            fprintf('- Invalid orbital elements\n');
        end
    end
end