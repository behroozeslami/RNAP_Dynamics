function [] = fit_pauses_Singletraces(row)

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

Pause_threshold =  37.88;

bins_per_decade = 7;

Ti = 0.01;

number_of_decades = 6;

correction = 0;

t_strt = 0.1;

t_end = 10000;

fitCount_tot = 3;

pause_num = 3;

Dobcktrck = 0;

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

Fit_folder = [output_foldername slash 'Fit_Results_ST_BS' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

if ~(exist(Fit_folder,'dir')==7)
    
    mkdir(Fit_folder)
    
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

dat_length_mat = load(str_DWT_length, '-mat');

dat_length_mat = dat_length_mat.dat_length;

dat_length = dat_length_mat(row);

DWT_dat = load(str_DWT_dat, '-mat');

DWT_dat = DWT_dat.DWT_dat;

DWT_array_0 = DWT_dat(row,1:dat_length);

DWT_hist_dat = load(str_DWT_hist, '-mat');

DWT_hist_dat = DWT_hist_dat.DWT_hist_dat;

DWT_hist_exp = DWT_hist_dat(row,:);

DWT_hist_err = load(str_DWT_error, '-mat');

DWT_hist_err = DWT_hist_err.DWT_hist_err;

selected_trace_number = load(str_trace_index, '-mat');

selected_trace_number = selected_trace_number.selected_trace_number;

trace_id = selected_trace_number(row);

Err_exp = DWT_hist_err(row,:);

[bar_pos, bins]=make_bins(Ti,bins_per_decade,number_of_decades,t0,correction);

gamma_nums = [Nw, ones([1,pause_num - 1])];

lb = [20,   14    0.15    0.005    0.01    0.001  10^-4    0.2    0.31];

ub = [20    25    0.45    0.015    0.2    0.02    2*10^-3    0.8    0.36];

params_0 = 0.5*(ub+lb);

options = saoptimset('TolFun',10^-10, 'MaxFunEvals', 10^5);

str_fit = [Fit_folder slash 'Fit_Result' '_' 'Trace' '_' num2str(trace_id) '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

PL = figure();

ax = axes();

for fitCount=1:fitCount_tot
    
    DWT_array = randsample(DWT_array_0,dat_length,true);
    
    [DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);

    binwidth = bins(2:length(bins)) - bins(1:(length(bins)-1));

    Weights = dat_length * DWT_hist .* binwidth;

    bool_vec = Weights > 0;

    Weights = Weights(bool_vec);

    bincenter = bar_pos(bool_vec);

    MLfun= @(params) ML_func_CG_nolog(params, B, Weights, bincenter, gamma_nums, Nw, t_strt, t_end, Dobcktrck);

    [fit,fval] = simulannealbnd(MLfun, params_0, lb, ub, options);
    
    k_array = fit(1:pause_num+1);

    q_array = fit(pause_num+2:2*pause_num+1);

    sigma_ln = fit(2*pause_num+2);

    t_bar = fit(2*pause_num+3);

    T_array=10.^(-1:0.01:3);

    y = Calculate_Model_v2(T_array, q_array, k_array, gamma_nums, sigma_ln, t_bar, Nw, B, t_strt, t_end, Dobcktrck);

    fit_dat = fit;
    
    save([str_fit, '.txt'], 'fit_dat', '-ascii', '-append')
                                                                                                                                                                                                                                        
    loglog(T_array, y, 'g', 'Linewidth', 0.5);
    
    hold on
    
end

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

close();

%exit;





