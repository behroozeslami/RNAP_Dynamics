function [fit_dat_All] = Analyze_fit_pauses_Singletraces(row)

if ischar(row)
    
    row = round(str2double((row)));
    
end

f0 = 12;

f = 7.5;

B = exp(-f/f0);

W = 12;

Nw = 10;

cam_freq = 25;

t0 = 1/cam_freq;

Pause_threshold =  4.84;

bins_per_decade = 7;

Ti = 0.01;

number_of_decades = 6;

correction = 0;

t_strt = 0.1;

t_end = 10000;

pause_num = 3;

Dobcktrck = 0;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'Export_AF_7.5pN_1mM_NTP_I';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

Fit_folder = [output_foldername slash 'Fit_Results_ST_BS' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

BestFit_folder = [output_foldername slash 'BestFit_ST_BS' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

if ~(exist(BestFit_folder,'dir')==7)
    
    mkdir(BestFit_folder)
    
end

DWT_folder = [output_foldername slash 'selected_DWT_pdf_Corrected' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_DWT = [DWT_folder slash 'DWT_pdf_Corrected' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

DWT_hist_dat = load([str_DWT, '.txt'], '-ascii');

bins_exp = transpose(DWT_hist_dat(:,1));

bar_pos_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,2));

DWT_hist_exp_comb = transpose(DWT_hist_dat(1:length(bins_exp)-1,3));

Err_exp_comb = transpose(DWT_hist_dat(1:length(bins_exp)-1,4));

str_DWT_length = [output_foldername slash 'DWT_data_length_Corrected' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's' '.mat'];

str_DWT_dat = [output_foldername slash 'DWT_data_Corrected' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's' '.mat'];

str_DWT_hist = [output_foldername slash 'DWT_pdf_Corrected' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's' '.mat'];

str_DWT_error = [output_foldername slash 'DWT_pdf_error_Corrected' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's' '.mat'];

str_trace_index = [output_foldername slash 'Trace_indices' '.mat'];

DWT_hist_dat = load(str_DWT_hist, '-mat');

DWT_hist_dat = DWT_hist_dat.DWT_hist_dat;

DWT_hist_exp = DWT_hist_dat(row,:);

DWT_hist_err = load(str_DWT_error, '-mat');

DWT_hist_err = DWT_hist_err.DWT_hist_err;

selected_trace_number = load(str_trace_index, '-mat');

selected_trace_number = selected_trace_number.selected_trace_number;

trace_id = selected_trace_number(row);

Err_exp = DWT_hist_err(row,:);

[bar_pos, bins] = make_bins(Ti,bins_per_decade,number_of_decades,t0,correction);

gamma_nums = [Nw, ones([1,pause_num - 1])];

str_fit = [BestFit_folder slash 'BestFit_Result' '_' 'Trace' '_' num2str(trace_id) '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_fit_All = [Fit_folder slash 'Fit_Result' '_' 'Trace' '_' num2str(trace_id) '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

fit_dat_All = load([str_fit_All, '.txt']', '-ascii');

bool_vec = (fit_dat_All(2,:) <2.6) & (fit_dat_All(4,:) > 5*10^-3) & (fit_dat_All(7,:) > 10^-4);

fit_dat_All = fit_dat_All(bool_vec,:);

if ~isempty(fit_dat_All)
    
    fit_dat = mean(fit_dat_All);

    k_array = fit_dat(1:pause_num+1);

    q_array = fit_dat(pause_num+2:2*pause_num+1);

    sigma_ln = fit_dat(2*pause_num+2);

    t_bar = fit_dat(2*pause_num+3);

    T_array=10.^(-1:0.01:3);

    y = Calculate_Model_v2(T_array, q_array, k_array, gamma_nums, sigma_ln, t_bar, Nw, B, t_strt, t_end, Dobcktrck);

    q_array(pause_num) = 0;

    y0 = Calculate_Model_v2(T_array, q_array, k_array, gamma_nums, sigma_ln, t_bar, Nw, B, t_strt, t_end, Dobcktrck);

    PL = figure();

    ax = axes();

    loglog(T_array, y, 'g', 'Linewidth', 2.5);

    hold on

    loglog(T_array, y0, 'b', 'Linewidth', 2);

    errorbar(bar_pos_exp ,DWT_hist_exp_comb, Err_exp_comb,'kO','MarkerSize',6, 'MarkerFaceColor', 'k');

    errorbar(bar_pos_exp ,DWT_hist_exp, Err_exp,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

    hold off

    ylim([10^-8, 5])

    str_title = ['Trace' num2str(trace_id) ' ' 'Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw)];

    title(str_title)

    xlabel('Dwell Time (s)')

    ylabel('PDF')

    set(ax,'XScale','log','YScale','log')

    saveas(PL, [str_fit, '.jpg'], 'jpg')

    saveas(PL, [str_fit, '.fig'], 'fig')

    save([str_fit, '.txt'], 'fit_dat', '-ascii')

    close All
    
end







