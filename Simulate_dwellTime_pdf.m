function [DWT_hist] = Simulate_dwellTime_pdf(k_array,q_array,Dobcktrck,t0,f,Trace_time,Trace_num,W,k,Nw,sigma,bins)

DWT_array = [];

for i=1:Trace_num
    
    disp(i)
    
    [t_array, x_array] = Multistate_Gillespie(k_array,q_array,Dobcktrck,t0,f,Trace_time);
    
    x_noise = gen_noise(sigma, t0, Trace_time);
    
    x_array = x_array + x_noise;
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);
    
    [Tc_array] = find_crossing_time(t_array, xs_array,Nw);
    
    [DWT_array_now]=find_Dwell_time(Tc_array);
        
    DWT_array = [DWT_array, DWT_array_now];
    
end

[DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);