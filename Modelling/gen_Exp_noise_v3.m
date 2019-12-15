function x_noise = gen_Exp_noise_v3(all_Pauses,points)

x_noise = zeros([1,points]);

Locs = 1:length(all_Pauses)-points+1;

index = 0;

while index < points
    
    Loc = randsample(Locs,1);
    
    Pause = all_Pauses(Loc:Loc+points-1);
    
    new_index = index + length(Pause);
    
    new_index = min(new_index, points);
    
    x_noise(index+1:new_index) = Pause(1:new_index-index);
    
    index = new_index; 
end

