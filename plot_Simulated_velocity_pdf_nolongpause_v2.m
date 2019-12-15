
k_array = [24.8388    1.6652    0.2559    0.0570    0.0041];

q_array = [0.0741    0.0062    0.0016    0.0004];

Dobcktrck = 0;

cam_freq = 25;

t0 = 1/cam_freq;

f0 = 12;

f = 7.5;

B = exp(-f/f0);

Trace_time = 10^4;

Trace_num = 20;

W = 12;

k = 0;

Nw = 10;

Ts = 2*W+1;

T = Ts;

Pause_threshold_V = 100.64;

binwidth = 1;

Vmax = 60;

Vmin = -40;

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'Export_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

noise_V_pdf_folder = [output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise_V_dat = [noise_V_pdf_folder slash 'NoiseVelocity' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

noise_velocities_sample = load([str_noise_V_dat '.txt'],'-ascii');

[centers, simul_velocity_pdf] = Simulate_velocity_pdf_nolongpause_v2(k_array,q_array,Dobcktrck,t0,B,Trace_time,Trace_num,W,k,Nw,T,Pause_threshold_V,noise_velocities_sample,binwidth,Vmin,Vmax);

figure()

PL = plot(centers, simul_velocity_pdf, 'c-', 'LineWidth', 2)

legend(['Simulation (k_{el} = ' num2str(k_array(1)) ' bp/s)'], 'Location', 'South')

xlabel('Velocity (bp/s)')

ylabel('PDF')