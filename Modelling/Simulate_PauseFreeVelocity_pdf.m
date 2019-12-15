function [V_hist, centers] = Simulate_PauseFreeVelocity_pdf(Model_ID, Params, t0, Trace_time, V_num, W, T, Nw, Pause_threshold_V, noise, noise_type, Vmin, Vmax, binwidth)
    
if strcmp(noise_type, 'Gaussian')
    
    Gen_noise_Func = @(points) gen_noise(noise, points);
    
end
    
if strcmp(noise_type, 'Exp')
        
    Gen_noise_Func = @(points) gen_Exp_noise_v3(noise, points);
        
end

if Model_ID == 1
    
    Gillespie_Func = @(Trace_time) Gillespie_MicroModel_1(Params, t0, Trace_time);
end
 
Dists = [];

sample_size = 0;

while sample_size < V_num
    
    %disp(sample_size)
    
    [t_array, x_array, s_array] = feval(Gillespie_Func,Trace_time);
    
    points = length(x_array);
    
    [Max, Ind] = max(x_array);
    
    points = min([Ind, points]);
    
    t_array =  t_array(1:points);
        
    x_array =  x_array(1:points);
    
    if length(t_array) > 2*W+1
        
        x_noise = feval(Gen_noise_Func, points);
    
        x_array = x_array + x_noise;
        
        [xs_array]=Smooth_Trace_SG(x_array,W, 0);
        
        len = length(xs_array);
    
        Dist_array = transpose(xs_array(T+1:len) - xs_array(1:len-T));

        [Tc_array] = find_crossing_time(t_array, xs_array, Nw);
    
        [DWT_array]=find_Dwell_time(Tc_array);
        
        DWT_indices = 1:length(DWT_array);
    
        pause_start_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V));
    
        pause_end_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V)+1);
    
        pause_start_indices = round(pause_start_times/t0 +1);
    
        pause_end_indices = round(pause_end_times/t0 +1);
    
        pause_start_indices = pause_start_indices-(T+1);
    
        pause_start_indices(pause_start_indices < 1) = 1;
    
        bool = true([1, len]);
        
        for j=1:length(pause_start_indices)
        
            bool(pause_start_indices(j):pause_end_indices(j)) = false(1);
        
        end
        
        bool = bool(1:length(Dist_array));
        
        Dist_array = Dist_array(bool);
        
        sample_size = sample_size+length(Dist_array);
        
        Dists = [Dists, Dist_array'];
        
    end
    
end

velocities = Dists/(T*t0);

centers = Vmin:binwidth:Vmax;

V_hist = hist(velocities, centers);

V_hist = V_hist/(binwidth*sum(V_hist));

centers = centers(2:length(centers)-1);

V_hist = V_hist(2:length(V_hist)-1);

