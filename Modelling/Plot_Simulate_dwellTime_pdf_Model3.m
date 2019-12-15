f = -12.5;

f0 = 12;

delta = 0.7;

W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 4;
     
T0 = (2*W+1);

Pause_threshold_V = 100.64;
%Pause_threshold_V =  37.88;

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0 = 1/cam_freq;

Model_ID = 3

str_Param = 'Parameters/Micro_Model_3_NewSel_Parameters_NewSel_AF_7.5pN_1mM_NTP.xlsx';


Trace_time = 3*10^3;

DWT_num = 10^5;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

%mainfoldername = 'Select_AF_7.5pN_1mM_NTP';
%mainfoldername = 'EXPORT_AF_7.5pN_0.5mM_NTP';
mainfoldername = 'Export_AF_7.5pN_1mM_NTP';
%mainfoldername = 'NEW_AF_7.5pN_10uM_NTP';
%mainfoldername = 'NEW_AF_7.5pN_1mM_ITP_100uM_NTP';
%mainfoldername = 'EXPORT_OF_12.5pN_1mM_NTP';
%mainfoldername = 'EXPORT_OF_5pN_1mM_NTP';
%mainfoldername = 'EXPORT_OF_7.5pN_1mM_NTP_NoB2';
%mainfoldername = 'NEW_AF_7.5pN_1mM_ITP_100uM_NTP';
%mainfoldername = 'EXPORT_OF_7.5pN_1mM_NTP_NoB2';
%mainfoldername = 'OF4_9pN_1mM_NTP_GreB_M&H';
%mainfoldername = 'NewSel_AF_7.5pN_1mM_NTP';


subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

DWT_folder = ['..' slash output_foldername slash 'selected_DWT_pdf_Combined' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

str_DWT = [DWT_folder slash 'DWT_pdf_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

str_DWT_dat = [DWT_folder slash 'DWT_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

DWT_hist_dat = load([str_DWT, '.txt'], '-ascii');

bins_exp = transpose(DWT_hist_dat(:,1));

bar_pos_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,2));

DWT_hist_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,3));

Err = transpose(DWT_hist_dat(1:length(bins_exp)-1,4));

noise_V_pdf_folder = ['..' slash output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T0) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise = [noise_V_pdf_folder slash 'NoiseSample_small' '_' 'Ts='  num2str(2*W+1) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '.mat'];

noise_dat = load(str_noise, '-mat');
    
noise = transpose(noise_dat.noise_dat_small);

noise_type = 'Exp';

[Param_dat, txt] = xlsread(str_Param);

Params = Param_dat(1,:);

%Params(1) = 17;

%Params(2) = 0.26;


force_indices = [6];

q = Params(force_indices(1));

Params(force_indices(1)) = q*exp((7.5-f)*(1-delta)/f0);

disp(Params)

Simul_DWT_hist = Simulate_dwellTime_pdf(Model_ID, Params, t0, Trace_time, DWT_num, W, Nw, noise, noise_type, bins_exp);

figure()

ax = axes();

PL1 = errorbar(bar_pos_exp,DWT_hist_exp,Err,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

hold on

loglog(bar_pos_exp, Simul_DWT_hist, 'g-', 'LineWidth', 2)

hold off

str_title = ['Ts = ' ' ' num2str((2*W+1)*t0) 's' ' ' 'Nw = ' num2str(Nw)];

title(str_title)

xlabel('Dwell Time (s)')

ylabel('PDF')

set(ax,'XScale','log','YScale','log')


