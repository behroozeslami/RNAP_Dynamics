% Update: I save the error for the parameters of the micro-model
% parameters. Also I correct for the kel and P (if the output file exists)

f0 = 12;

f = 7.5;

B = exp(-f/f0);

W = 12;

k = 0;
  
Nw = 10;

cam_freq = 25;

t0 = 1/cam_freq;

%Pause_threshold =  37.88;
Pause_threshold =  100.64;

bins_per_decade = 7;

Ti = 0.01;

number_of_decades = 6;

correction = 0;

t_strt = 0.1;

t_end = 10000;

fitCount_tot_1 = 10;

fitCount_tot_2 = 50;

pause_num = 3;

Dobcktrck = 0;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'NewSel_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

Fit_folder = [output_foldername slash 'Fit_Results_BS' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

if ~(exist(Fit_folder,'dir')==7)
    
    mkdir(Fit_folder)
    
end


DWT_folder = [output_foldername slash 'selected_DWT_pdf_Corrected' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_DWT = [DWT_folder slash 'DWT_pdf_Corrected' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_DWT_dat = [DWT_folder slash 'DWT_Corrected' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

DWT_hist_dat = load([str_DWT, '.txt'], '-ascii');

bins_exp = transpose(DWT_hist_dat(:,1));

bar_pos_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,2));

DWT_hist_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,3));

Err_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,4));

DWT_array_0 = transpose(load([str_DWT_dat, '.txt'], '-ascii'));

DWT_array_0 = DWT_array_0(DWT_array_0 > t_strt & DWT_array_0 < t_end);

dat_length = length(DWT_array_0);

[bar_pos, bins]=make_bins(Ti,bins_per_decade,number_of_decades,t0,correction);

gamma_nums = [Nw, ones([1,pause_num - 1])];
%gamma_nums = Nw*ones([1,pause_num])

Blacklist = load([Fit_folder slash 'Blacklist_orig.txt']);

indices2keep = true(1, fitCount_tot_1);

indices2keep(Blacklist) = false;

fitCounts = 1:fitCount_tot_1;

fitCounts = fitCounts(indices2keep);

fitNum = length(fitCounts);

Min_ML = Inf;

for j=1:fitNum
    
    fitCount = fitCounts(j);
    
    str_fit = [Fit_folder slash 'Fit_Result_OrigDat' '_' num2str(fitCount) '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];
    
    fit_dat = load([str_fit, '.txt'], '-ascii');
    
    model_params = transpose(fit_dat(:,1));
    
    ML = fit_dat(1,3);
    
    if ML < Min_ML
        
        Min_ML = ML;
        
        best_model_params = model_params;
        
        best_fitCount = fitCount;
        
    end
end

Blacklist = load([Fit_folder slash 'Blacklist.txt']);

indices2keep = true(1, fitCount_tot_2);

indices2keep(Blacklist) = false;

fitCounts = 1:fitCount_tot_2;

fitCounts = fitCounts(indices2keep);

fitNum = length(fitCounts);

sum_SQ_delta = zeros([1, 2*pause_num+3]);

All_model_params = zeros([fitNum,  2*pause_num+3]);

for j=1:fitNum
    
    fitCount = fitCounts(j);
    
    str_fit = [Fit_folder slash 'Fit_Result' '_' num2str(fitCount) '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];
    
    fit_dat = load([str_fit, '.txt'], '-ascii');
    
    model_params = transpose(fit_dat(:,1));
    
    All_model_params(j,:) = model_params';
    
    sum_SQ_delta = sum_SQ_delta + (model_params - best_model_params).^2;
end

Err_model_params = sqrt(sum_SQ_delta/fitNum);

Low_Bond = quantile(All_model_params,0.136,1);

Up_Bond = quantile(All_model_params,1-0.136,1);

mean_model_params = mean(All_model_params);

Conf_int_dat = [best_model_params; mean_model_params; Low_Bond; Up_Bond];

k_array = best_model_params(1:(pause_num+1));

q_array = best_model_params((pause_num+2):(2*pause_num+1));

t_bar = best_model_params(2*pause_num+2);

sigma_ln = best_model_params(2*pause_num+3);

Err_k_array = Err_model_params(1:pause_num+1);

Err_q_array = Err_model_params((pause_num+2):(2*pause_num+1));

len = pause_num+1;

best_model_param_dat = NaN(5, len);

best_model_param_dat(1,:) = k_array;

best_model_param_dat(2,:) = Err_k_array;

best_model_param_dat(3,1:(len-1)) = q_array;

best_model_param_dat(4,1:(len-1)) = Err_q_array;

best_model_param_dat(5,1) = t_bar;

best_model_param_dat(5,2) = sigma_ln;

T_array=10.^(-1:0.01:3.5);

y = Calculate_Model_v2(T_array, q_array, k_array, gamma_nums, sigma_ln, t_bar, Nw, B, t_strt, t_end, Dobcktrck);

figure()

ax = axes();

PL = loglog(T_array, y, 'g', 'Linewidth', 2);

set(ax,'XScale','log','YScale','log')

hold on

errorbar(bar_pos_exp ,DWT_hist_exp, Err_exp,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');
    
hold off

ylim([10^-8, 5])

str_title = ['Ts =' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

xlabel('Dwell Time (s)')

ylabel('PDF')

str_best_fit = [Fit_folder slash 'Best_Pause_Fit' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

save([str_best_fit, '.txt'], 'best_model_param_dat', '-ascii')

saveas(PL, [str_best_fit, '.jpg'], 'jpg')

saveas(PL, [str_best_fit, '.fig'], 'fig')

str_best_fit_confint = [Fit_folder slash 'Best_Pause_Fit_ConfInt' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

save([str_best_fit_confint, '.txt'], 'Conf_int_dat', '-ascii')

dat = Conf_int_dat(:, 1:7);

param_foldername = ['Modelling' slash 'Parameters'];

if ~(exist(param_foldername,'dir')==7)
    
    mkdir(param_foldername)
    
end

output_file_str = [param_foldername slash 'General_Model_NewSel_Parameters' '_' output_foldername '.xlsx'];

col_header = {'kel', 'k1', 'k2', 'k3', 'P1', 'P2', 'P3'}; 

row_header(1:4,1) = {'Best', 'Mean', 'LB', 'UB'};

xlswrite(output_file_str,dat,'Sheet1','B2');

xlswrite(output_file_str,col_header,'Sheet1','B1');

xlswrite(output_file_str,row_header,'Sheet1','A2'); 

All_micromodel_params = cell2mat(arrayfun(@(row) Micromodel_1_params(All_model_params(row,1:pause_num+1), All_model_params(row,pause_num+2:2*pause_num+1))', 1:length(All_model_params), 'UniformOutput', false ))';

dat_micro = [mean(All_micromodel_params); std(All_micromodel_params)];

file_to_correct_kel_P = './Modelling/Export_AF_7.5pN_1mM_NTP/Analysis_1mM_SingleTraces_kel-P/1mM_MV_2/Fit_outcome_1mM_MV_2_kel-P_Ts=1s_Nw=4_limits=0.05-10000.xlsx';

if exist(file_to_correct_kel_P, 'file') == 2
    
    correct_dat = xlsread(file_to_correct_kel_P);
    
    dat_micro(1,1:2) = correct_dat(1,:);
    
    dat_micro(2,1:2) = correct_dat(4,:);
    
end

output_file_str_2 = [param_foldername slash 'Micro_Model_1_NewSel_Parameters' '_' output_foldername '.xlsx'];

row_header_2(1:2,1) = {'Mean', 'std'};

col_header_2 = {'kel', 'P', 'ke1', 'ke2', 'ke3', 'kp2', 'kp3'}; 

row_header_2(1,1) = {'Mean'};

xlswrite(output_file_str_2,dat_micro,'Sheet1','B2');

xlswrite(output_file_str_2,col_header_2,'Sheet1','B1');

xlswrite(output_file_str_2,row_header_2,'Sheet1','A2');  



