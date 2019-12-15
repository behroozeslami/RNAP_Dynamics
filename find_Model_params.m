f0 = 12;

f = 7.5;

B = exp(-f/f0);

W = 12;

k = 0;
  
Nw = 10;

cam_freq = 25;

t0 = 1/cam_freq;

Pause_threshold =  4.84;

t_strt = 0.1;

t_end = 10000;

fitCount_tot_2 = 50;

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

Fit_folder = [output_foldername slash 'Fit_Results_BS' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

Blacklist = load([Fit_folder slash 'Blacklist.txt']);

indices2keep = true(1, fitCount_tot_2);

indices2keep(Blacklist) = false;

fitCounts = 1:fitCount_tot_2;

fitCounts = fitCounts(indices2keep);

fitNum = length(fitCounts);

All_model_params = zeros([fitNum,  2*pause_num+3]);

for j=1:fitNum
    
    fitCount = fitCounts(j);
    
    str_fit = [Fit_folder slash 'Fit_Result' '_' num2str(fitCount) '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];
    
    fit_dat = load([str_fit, '.txt'], '-ascii');
    
    model_params = transpose(fit_dat(:,1));
    
    All_model_params(j,:) = model_params';
    
end

kel_dist = load([Fit_folder slash 'kel_Prob_grid_3Pauses_lowP.txt'], '-ascii');

P1_dist = load([Fit_folder slash 'P1_Prob_grid_3Pauses_lowP.txt'], '-ascii');

All_model_params(:, 5) = transpose(arrayfun(@(row) genDist(P1_dist(:,2), P1_dist(:,1)), 1:length(All_model_params)));

All_model_params(:, 1) = transpose(arrayfun(@(row) genDist(kel_dist(:,2), kel_dist(:,1)), 1:length(All_model_params)));

Low_Bond = quantile(All_model_params,0.136,1);

Up_Bond = quantile(All_model_params,1-0.136,1);

mean_model_params = mean(All_model_params);

Conf_int_dat = [mean_model_params; Low_Bond; Up_Bond];

dat = Conf_int_dat(:, 1:7);

param_foldername = ['Modelling' slash 'Parameters'];

if ~(exist(param_foldername,'dir')==7)
    
    mkdir(param_foldername)
    
end

output_file_str = [param_foldername slash 'General_Model_Parameters' '_' output_foldername '.xlsx'];

col_header = {'kel', 'k1', 'k2', 'k3', 'P1', 'P2', 'P3'}; 

row_header(1:3,1) = {'Mean', 'LB', 'UB'};

xlswrite(output_file_str,dat,'Sheet1','B2');

xlswrite(output_file_str,col_header,'Sheet1','B1');

xlswrite(output_file_str,row_header,'Sheet1','A2');  

Micro_params = transpose(cell2mat((arrayfun(@(row) Micromodel_1_params(All_model_params(row, 1:pause_num+1), All_model_params(row, pause_num+2:2*pause_num+1))', 1:length(All_model_params), 'UniformOutput', false))));

Low_Bond_Micro = quantile(Micro_params,0.136,1);

Up_Bond_Micro = quantile(Micro_params,1-0.136,1);

mean_model_params_Micro = mean(Micro_params);

Micro_dat = [mean_model_params_Micro; Low_Bond_Micro; Up_Bond_Micro];

Micro_col_header = {'kel', 'P', 'ke1', 'ke2', 'ke3', 'kp2', 'kp3'}; 

output_file_str2 = [param_foldername slash 'Micro_Model_1_Parameters' '_' output_foldername '.xlsx'];

xlswrite(output_file_str2,Micro_dat,'Sheet1','B2');

xlswrite(output_file_str2,Micro_col_header,'Sheet1','B1');

xlswrite(output_file_str2,row_header,'Sheet1','A2');  



