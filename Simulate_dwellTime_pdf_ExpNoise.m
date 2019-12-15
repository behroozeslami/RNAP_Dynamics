function [DWT_hist] = Simulate_dwellTime_pdf_ExpNoise(k_array,q_array,Dobcktrck,t0,B,Trace_time,Trace_num,W,k,Nw,noise_dat,bins)

DWT_array = [];

for i=1:Trace_num
    
    disp(i)
    
    [t_array, x_array] = Multistate_Gillespie(k_array,q_array,Dobcktrck,t0,B,Trace_time);
    
    [xs_array] = Smooth_Trace_SG(x_array,W, k);
    
    %xs_noise = gen_Exp_noise(noise_dat, t0, Trace_time);
    
    xs_noise = gen_Exp_noise_v2(noise_dat,W, k, t0, Trace_time);
    
    xs_array = xs_array + xs_noise;
    
    [Tc_array] = find_crossing_time(t_array, xs_array,Nw);
    
    [DWT_array_now]=find_Dwell_time(Tc_array);
        
    DWT_array = [DWT_array, DWT_array_now];
    
end

[DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);