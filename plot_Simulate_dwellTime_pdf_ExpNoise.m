cam_freq = 25;

t0 = 1/cam_freq;

f0 = 12;

f = 7.5;

B = exp(-f/f0);

Trace_time = 10^4;

Trace_num = 100;

W = 12;

Ts = (2*W+1);

T = Ts;

k = 0;

Nw = 10;

Pause_threshold =  37.88;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

Dobcktrck = 0;

t_strt = 0.1;

t_end = 10000;

mainfoldername = 'Select_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

Simul_DWT_folder = [output_foldername slash 'Simulated_DWT_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

if ~(exist(Simul_DWT_folder,'dir')==7)
    
    mkdir(Simul_DWT_folder)
    
end

str_Simul_DWT = [Simul_DWT_folder slash 'Simulated_DWT_pdf_2' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];
%{
Best_Fit_folder = [output_foldername slash 'Best_Fit' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_k_array = [Best_Fit_folder slash 'BestFit_k_array' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'T=' num2str(T) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_q_array = [Best_Fit_folder slash 'BestFit_q_array' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'T=' num2str(T) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

k_dat = load([str_k_array, '.txt'], '-ascii');

k_array = k_dat(1,:);

q_dat = load([str_q_array, '.txt'], '-ascii');
%}

k_array = [ 23.0144    1.8936    0.2643    0.0093];

q_array = [0.07    0.007    0.0007];

noise_V_pdf_folder = [output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise = [noise_V_pdf_folder slash 'NoiseSample' '_' 'Ts='  num2str(2*W+1) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

noise_dat = load([str_noise, '.txt'], '-ascii');

DWT_folder = [output_foldername slash 'selected_DWT_pdf_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

str_DWT = [DWT_folder slash 'DWT_pdf_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

DWT_hist_dat = load([str_DWT, '.txt'], '-ascii');

bins = transpose(DWT_hist_dat(:,1));

bar_pos = transpose(DWT_hist_dat(1:length(bins)-1,2));

DWT_hist_exp = transpose(DWT_hist_dat(1:length(bins)-1,3));

Err = transpose(DWT_hist_dat(1:length(bins)-1,4));

DWT_folder_correct = [output_foldername slash 'selected_DWT_pdf_Corrected' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_DWT_correct = [DWT_folder_correct slash 'DWT_pdf_Corrected' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

DWT_hist_dat_correct = load([str_DWT_correct, '.txt'], '-ascii');

DWT_hist_exp_correct = transpose(DWT_hist_dat_correct(1:length(bins)-1,3));

Err_correct = transpose(DWT_hist_dat_correct(1:length(bins)-1,4));

[Simul_DWT_hist] = Simulate_dwellTime_pdf_ExpNoise(k_array,q_array,Dobcktrck,t0,B,Trace_time,Trace_num,W,k,Nw,noise_dat,bins);

[Simul_DWT_hist_zero_noise] = Simulate_dwellTime_pdf(k_array,q_array,Dobcktrck,t0,B,Trace_time,Trace_num,W,k,Nw,0,bins);

figure()

ax = axes();

PL = errorbar(bar_pos,DWT_hist_exp,Err,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

hold on

errorbar(bar_pos,DWT_hist_exp_correct,Err_correct,'rO','MarkerSize',6);

loglog(bar_pos, Simul_DWT_hist, 'k-', 'LineWidth', 2)

loglog(bar_pos, Simul_DWT_hist_zero_noise, 'b--', 'LineWidth', 1.5)

hold off

set(ax,'XScale','log','YScale','log')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

legend('Experimental distribution (not corrected)', 'Experimental distribution (corrected)', 'Simulation (with experimental noise)', 'Simulation (filtered traces, zero noise)', 'Location','southWest')

xlabel('Dwell Time (s)')

ylabel('PDF')

saveas(PL, [str_Simul_DWT, '.jpg'], 'jpg')

saveas(PL, [str_Simul_DWT, '.fig'], 'fig')