 function [Fit_dat, orig_simul_velocity_pdf] = Calculate_Velocity_LogLikelihood_bootstrapping(Params, V_Mat, Vmin_fit, Vmax_fit, Model_ID, t0, Trace_time, V_num, W, T, Nw, Pause_threshold_V, noise, noise_type, Vmin, Vmax, binwidth)
 
Ts = 2*W+1;

Size = size(V_Mat);

nboot = Size(1);

dat_length = Size(2);

[simul_velocity_pdf, centers] = Simulate_PauseFreeVelocity_pdf(Model_ID, Params, t0, Trace_time, V_num, W, T, Nw, Pause_threshold_V, noise, noise_type, Vmin, Vmax, binwidth);

orig_simul_velocity_pdf = simul_velocity_pdf;

Norm = sum(simul_velocity_pdf((centers > Vmin_fit)&(centers < Vmax_fit))*binwidth); 
 
simul_velocity_pdf = simul_velocity_pdf/Norm;

log_simul_velocity_pdf = log(simul_velocity_pdf);
        
log_simul_velocity_pdf(abs(log_simul_velocity_pdf)==Inf) = 10^(-300);

centers = Vmin:binwidth:Vmax;

LogLike_Func = @(row) -(binwidth*dat_length/Ts) * dot(V_hist(V_Mat(row,:), centers, binwidth), log_simul_velocity_pdf);

Fit_dat = arrayfun(LogLike_Func, 1:nboot);

    
    
    