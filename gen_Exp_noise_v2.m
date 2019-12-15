function xs_noise = gen_Exp_noise_v2(noise_dat,W, k, t0, Trace_time)

Ts = 2*W+1;

max_index = floor(Trace_time/t0) + 1;

all_Pauses = transpose(noise_dat(:,1));

Pause_number_array = transpose(noise_dat(:,3));

Pause_num_tot = max(Pause_number_array);

Pause_nums = 1:Pause_num_tot;

xs_noise = [];

while length(xs_noise) < max_index
    
    Pause_num = randsample(Pause_nums,1, true);
    
    Pause = all_Pauses(Pause_number_array == Pause_num);
    
    len = length(Pause);
    
    if len > 2*Ts
        
        Pause_s = Smooth_Trace_SG(Pause,W, k);
        
        Pause_s = Pause_s(W:len-W);
        
        xs_noise = [xs_noise, Pause_s];
        
    end
    
end

xs_noise = xs_noise(1:max_index);

