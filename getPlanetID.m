function p_id = getPlanetID(segment_num)
    % Maps segment number to planet ID for E-V-E-V-E-M sequence
    sequence = [3,2,3,2,3,4];  % Earth(3)-Venus(2)-Earth-Venus-Earth-Mars(4)
    p_id = sequence(mod(segment_num-1,length(sequence))+1);
end