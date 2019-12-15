 function [Fit_dat, Simul_DWT_array, Proc_data] = Calculate_LogLikelihood_SingleTrace(Params, dat_length, DWT_Mat, t_strt, t_end, Model_ID, t0, Trace_time, DWT_num, W, Nw, noise, noise_type)
    
    bins_per_decade=7;

    Ti=0.01;

    number_of_decades=6;

    correction=0;

    [bar_pos, bins] = make_bins(Ti,bins_per_decade,number_of_decades,t0,correction);

    binwidth = bins(2:length(bins)) - bins(1:(length(bins)-1));
    
    select_bins = (bar_pos > t_strt) & (bar_pos < t_end);
    
    selected_binwidth = binwidth(select_bins);
            
    [Simul_DWT_hist, Simul_DWT_array, Proc_data]  = Simulate_dwellTime_pdf(Model_ID, Params, t0, Trace_time, DWT_num, W, Nw, noise, noise_type, bins);
            
    Norm = sum(Simul_DWT_hist(select_bins).* selected_binwidth); 
    
    Simul_DWT_hist = Simul_DWT_hist/Norm;
        
    log_simul_DWT_hist = log(Simul_DWT_hist);
        
    log_simul_DWT_hist(abs(log_simul_DWT_hist)==Inf) = 10^(-300);
        
    LogLike_Func = @(row) -dat_length(row) * dot((DwellTimeHist_v3(DWT_Mat(row,1:dat_length(row)), t0, bins).*binwidth), log_simul_DWT_hist);
    
    L = length(dat_length);

    Fit_dat = arrayfun(LogLike_Func, 1:L);
    
    
    
    