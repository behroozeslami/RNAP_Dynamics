function xs_noise = gen_Exp_noise(noise_dat, t0, Trace_time)

max_index = floor(Trace_time/t0) + 1;

all_Pauses = transpose(noise_dat(:,2));

Pause_number_array = transpose(noise_dat(:,3));

Pause_num_tot = max(Pause_number_array);

Pause_nums = 1:Pause_num_tot;

xs_noise = [];

while length(xs_noise) < max_index
    
    Pause_num = randsample(Pause_nums,1, true);
    
    Pause = all_Pauses(Pause_number_array == Pause_num);
    
    xs_noise = [xs_noise, Pause];
    
end

xs_noise = xs_noise(1:max_index);

