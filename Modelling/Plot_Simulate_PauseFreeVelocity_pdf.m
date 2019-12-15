W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 10;

Pause_threshold_V = 100.64;

cam_freq = 25;

t0 = 1/cam_freq;

Model_ID = 1;

str_Param = 'Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx';

Trace_time = 3*10^3;

V_num = 10^7;

binwidth = 1;

Vmax = 60;

Vmin = -40;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'Select_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

V_pdf_folder = ['..' slash output_foldername slash 'Velocity_pdf_NoLongPause_v2' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_V_pdf = [V_pdf_folder slash 'PauseFreeVelocity_PDF' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

V_hist_dat = load([str_V_pdf, '.txt'], '-ascii');

centers = V_hist_dat(:,1);

Exp_V_hist = V_hist_dat(:,2);

noise_V_pdf_folder = ['..' slash output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise = [noise_V_pdf_folder slash 'NoiseSample_small' '_' 'Ts='  num2str(2*W+1) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '.mat'];

noise_dat = load(str_noise, '-mat');
    
noise = transpose(noise_dat.noise_dat_small);

noise_type = 'Exp';

[Param_dat, txt] = xlsread(str_Param);

Params = Param_dat(1,:);

Params(1) = 21.64;

Params(2) = 0.062;

disp(Params)

Simul_V_hist = Simulate_PauseFreeVelocity_pdf(Model_ID, Params, t0, Trace_time, V_num, W, T, Nw, Pause_threshold_V, noise, noise_type, Vmin, Vmax, binwidth);

figure()

PL1 = plot(centers,Exp_V_hist,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

hold on

loglog(centers, Simul_V_hist, 'b-', 'LineWidth', 2)

hold off

str_title = ['Ts = ' ' ' num2str((2*W+1)*t0) 's' ' ' 'T = ' num2str(T*t0)];

title(str_title)

xlabel('Velocity (bp/s)')

ylabel('PDF')

