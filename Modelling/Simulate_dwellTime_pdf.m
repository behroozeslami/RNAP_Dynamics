function [DWT_hist, DWT_array, Proc_data] = Simulate_dwellTime_pdf(Model_ID, Params, t0, Trace_time, DWT_num, W, Nw, noise, noise_type, bins)

DWT_array = zeros([1, DWT_num]);
    
if strcmp(noise_type, 'Gaussian')
    
    Gen_noise_Func = @(points) gen_noise(noise, points);
    
end
    
if strcmp(noise_type, 'Exp')
        
    Gen_noise_Func = @(points) gen_Exp_noise_v3(noise, points);
        
end

if Model_ID == 0
    
    Gillespie_Func = @(Trace_time) Gillespie_GeneralModel_3Pauses(Params, t0, Trace_time);
end
 
if Model_ID == 1
    
    Gillespie_Func = @(Trace_time) Gillespie_MicroModel_1(Params, t0, Trace_time);
end

if Model_ID == 2
    
    Gillespie_Func = @(Trace_time) Gillespie_MicroModel_2(Params, t0, Trace_time);
end

if Model_ID == 3
    
    Gillespie_Func = @(Trace_time) Gillespie_MicroModel_3(Params, t0, Trace_time);
end

index = 0;

Trace_num = 0;

sum_Proc = 0;

sum_Proc2 = 0;

while index < DWT_num
    
    %disp(index)
    
    [t_array, x_array, s_array] = feval(Gillespie_Func,Trace_time);
    
    points = length(x_array);
    
    [Max, Ind] = max(x_array);
    
    points = min([Ind, points]);
    
    t_array =  t_array(1:points);
        
    x_array =  x_array(1:points);
    
    if length(t_array) > 2*W+1
        
        Trace_num = Trace_num + 1;
        
        sum_Proc = sum_Proc + max(x_array);
        
        sum_Proc2 = sum_Proc2 + max(x_array)^2;
        
        x_noise = feval(Gen_noise_Func, points);
    
        x_array = x_array + x_noise;
        
        [xs_array]=Smooth_Trace_SG(x_array,W, 0);

        [Tc_array] = find_crossing_time(t_array, xs_array, Nw);
    
        [DWT_array_now]=find_Dwell_time(Tc_array);
    
        new_index = index + length(DWT_array_now);
    
        new_index = min(new_index, DWT_num);
        
        DWT_array(index+1:new_index) = DWT_array_now(1:new_index-index);
    
        index = new_index;
        
    end
    
end

[DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);

Proc = sum_Proc/Trace_num;

sigma_proc = sqrt((sum_Proc2/Trace_num-Proc^2))/sqrt(Trace_time);

Proc_data = [Proc, sigma_proc];