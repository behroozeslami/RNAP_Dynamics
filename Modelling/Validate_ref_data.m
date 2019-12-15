f = -5;

f0 = 12;

delta = 0.7;

W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 20;

orig_Nw = 10;
     
T0 = (2*W+1);

Pause_threshold_V =  37.88;

cam_freq = 25;

t0 = 1/cam_freq;

Model_ID = 1;

str_Param = 'Parameters/Micro_Model_1_Parameters_Select_AF_7.5pN_1mM_NTP.xlsx';

Trace_time = 3*10^3;

DWT_num = 10^5;

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

DWT_folder = ['..' slash output_foldername slash 'selected_DWT_pdf_Combined' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

str_DWT = [DWT_folder slash 'DWT_pdf_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

str_DWT_dat = [DWT_folder slash 'DWT_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

DWT_hist_dat = load([str_DWT, '.txt'], '-ascii');

bins_exp = transpose(DWT_hist_dat(:,1));

bar_pos_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,2));

DWT_hist_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,3));

Err = transpose(DWT_hist_dat(1:length(bins_exp)-1,4));

noise_V_pdf_folder = ['..' slash output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T0) '_' 'Nw=' num2str(orig_Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise = [noise_V_pdf_folder slash 'NoiseSample_small' '_' 'Ts='  num2str(2*W+1) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(orig_Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '.mat'];

noise_dat = load(str_noise, '-mat');
    
noise = transpose(noise_dat.noise_dat_small);

noise_type = 'Exp';

[Param_dat, txt] = xlsread(str_Param);

Params = Param_dat(1,:);

disp(Params)

Simul_DWT_hist = Simulate_dwellTime_pdf(Model_ID, Params, t0, Trace_time, DWT_num, W, Nw, noise, noise_type, bins_exp);

PL = figure();

ax = axes();

errorbar(bar_pos_exp,DWT_hist_exp,Err,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

hold on

loglog(bar_pos_exp, Simul_DWT_hist, 'k-', 'LineWidth', 2)

hold off

str_title = ['Ts = ' ' ' num2str((2*W+1)*t0) 's' ' ' 'Nw = ' num2str(Nw)];

title(str_title)

xlabel('Dwell Time (s)')

ylabel('PDF')

set(ax,'XScale','log','YScale','log')

header = {'Bin edges', 'Bin centers', 'Exp DWT pdf', 'Errors', 'Best fit (combined)'};

DWT_Simul_dat = [DWT_hist_dat, [Simul_DWT_hist, NaN(1,1)]'];

DWT_savefolder = [output_foldername slash 'Validation'];

if ~(exist(DWT_savefolder,'dir')==7)
    
    mkdir(DWT_savefolder)
    
end

str_DWT_save = [DWT_savefolder slash 'Validation_RefData_DWTpdf' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw)];
    
xlswrite([str_DWT_save,'.xlsx'],DWT_Simul_dat,'Sheet1','A2');

xlswrite([str_DWT_save,'.xlsx'],header,'Sheet1','A1');

saveas(PL,[str_DWT_save,'.fig'], 'fig')

saveas(PL,[str_DWT_save,'.jpg'], 'jpg')



