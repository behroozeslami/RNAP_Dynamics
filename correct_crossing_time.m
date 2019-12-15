function new_Tc_array = correct_crossing_time(Tc_array, Pause_threshold, t0)

%No_pause_prob = 1-Cum_pause_prob;

DWT_indices = 1:(length(Tc_array)-1);

indices2keep = true(1, length(Tc_array));

[DWT_array] = find_Dwell_time(Tc_array);

Pause_indices = DWT_indices((DWT_array- t0/2) >= Pause_threshold);

for j=2: length(Pause_indices)
    
    index = Pause_indices(j);
    
    prev_index = index - 1;
    
    if ((DWT_array(prev_index)-t0/2) >= Pause_threshold)
        
         indices2keep(index) = false;
        
    end
        
end
    
Tc_array = Tc_array(indices2keep);

new_Tc_array = Tc_array;