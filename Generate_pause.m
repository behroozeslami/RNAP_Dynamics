function Pause_tot = Generate_pause(Pause_sample, Pause_num)

Pauses = randsample(Pause_sample, Pause_num, 1);

Pause_tot = sum(Pauses);