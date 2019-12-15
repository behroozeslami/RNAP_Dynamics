function [centers, velocity_pdf] = Simulate_velocity_pdf_nolongpause_v2(k_array,q_array,do_bt,t0,B,Trace_time,Trace_num,W,k,Nw,T,Pause_threshold,noise_velocities_sample,binwidth,Vmin,Vmax)

Dists = [];

for i=1:Trace_num
    
    [t_array, x_array] = Multistate_Gillespie(k_array,q_array,do_bt,t0,B,Trace_time);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);
    
    [Tc_array] = find_crossing_time(t_array, xs_array,Nw);
    
    [DWT_array]=find_Dwell_time(Tc_array);
    
    DWT_indices = 1:length(DWT_array);
    
    pause_start_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold));
    
    pause_end_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold)+1);
     
    pause_start_indices = round(pause_start_times/t0 +1);
    
    pause_end_indices = round(pause_end_times/t0 +1);
    
    pause_bool_vec = ones([1, length(t_array)]);
    
    for n=1:length(pause_start_indices)
        
        pause_start = pause_start_indices(n);
        
        pause_end = pause_end_indices(n);
        
        pause_bool_vec(pause_start:pause_end) = 0;
    end
    
    delta_pause_bool_vec = pause_bool_vec(2:length(pause_bool_vec)) - pause_bool_vec(1:length(pause_bool_vec)-1);
    
    indices = 1:length(t_array);
    
    chunk_start = indices(delta_pause_bool_vec == 1);
    
    chunk_end = indices(delta_pause_bool_vec == -1);
    
    if (length(chunk_start) > 0) && (length(chunk_end) > 0)  
        
        if chunk_start(1) > chunk_end(1)
        
            chunk_start =  [1, chunk_start];
        end
        
        if chunk_start(length(chunk_start)) > chunk_end(length(chunk_end))
        
            chunk_end =  [chunk_end, indices(length(indices))];
        end
        
    else
        
        chunk_start = 1;
        
        chunk_end = length(t_array);
        
    end
    
    for j=1:length(chunk_start)
        
        x_chunk = xs_array(chunk_start(j):chunk_end(j));
        
        len = length(x_chunk);
    
        Dist_array = x_chunk(T+1:len) - x_chunk(1:len-T);
        
        Dists = [Dists, Dist_array];
        
    end
    
end

velocities = Dists/(T*t0);

len = length(velocities);

noise_velocities = randsample(noise_velocities_sample,len,true)';

velocities = velocities + noise_velocities;

centers = Vmin:binwidth:Vmax;

velocity_pdf = V_hist(velocities, centers, binwidth);

centers = centers(2:length(centers)-1);

velocity_pdf = velocity_pdf(2:length(velocity_pdf)-1);