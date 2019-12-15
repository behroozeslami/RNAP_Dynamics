function Analyze_1mM_kp2_SingleTraces(data_name, subfolder_name, selected_trace_indices) 

mainfoldername = 'Select_AF_7.5pN_1mM_NTP';

grid_name = 'ML_grid_1mM';

save_folder = 'Analysis_1mM_SingleTraces_kp2';

fit_index = 6;

param_num = 1;

t_strt = 0.05;

t_end = 30;

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 10;

alpha = 0.136; % 1-sigma confidence interval

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

save_subfolder = [output_foldername '/' save_folder '/' subfolder_name];

if ~(exist(save_subfolder,'dir')==7)
    
    mkdir(save_subfolder)
    
end

trace_indices = load([output_foldername slash 'Trace_indices.mat']);

trace_indices = trace_indices.selected_trace_number;

subset_cols =ismember(trace_indices, selected_trace_indices);

ndat = length(trace_indices);

[Param_dat, txt] = xlsread('Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Param_name = txt(1,fit_index+1);

folder_grid = [output_foldername slash save_folder];

str_grid = [folder_grid slash grid_name '_' Param_name{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

grid_dat = load(str_grid,'-mat');

grid = grid_dat.grid;

str_datlength = [output_foldername '/' 'DWT_data_length' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

DL = load(str_datlength, '-mat');

dat_length = DL.dat_length;

temp_grid = grid(:,param_num+1:ndat+param_num);

temp_grid = temp_grid(:,subset_cols);

grid = [grid(:,1:param_num), temp_grid];

dat_length = dat_length(subset_cols);

trace_indices = trace_indices(subset_cols);

ndat = length(trace_indices);

[Min_LL, best_indices] = min(grid(:,param_num+1:ndat+param_num));

gSize = size(grid);

points_num = gSize(1);

Prob_grid = zeros([points_num, ndat]);

for n=1:ndat
    
    ind = best_indices(n);
    
    Prob_grid(ind,n) = 1;
    
end

dat_length_normal = dat_length/sum(dat_length);

Ave_Prob_grid = Prob_grid*transpose(dat_length_normal);

kp2_bins = unique(grid(:,1));

delta_kp2 = kp2_bins(2)-kp2_bins(1);

[Best_kp2, Mean_kp2, LB_kp2, UB_kp2, std_kp2] = Conf_int(kp2_bins,Ave_Prob_grid,alpha);

Fit_dat = [Mean_kp2, LB_kp2, UB_kp2, std_kp2/sqrt(length(selected_trace_indices))]';

str_title = strrep(data_name, '_', ' - ');

PL1 = figure(1);

bar(kp2_bins,Ave_Prob_grid,'r')

hold on

line([LB_kp2-0.5*delta_kp2,LB_kp2-0.5*delta_kp2],ylim,'LineWidth', 2, 'Color', 'b')

line([UB_kp2+0.5*delta_kp2,UB_kp2+0.5*delta_kp2],ylim,'LineWidth', 2, 'Color', 'b')

line([Mean_kp2,Mean_kp2],ylim,'LineWidth', 2, 'Color', 'g')

hold off

xlabel('Optimal kp2 (bp/s)')

ylabel('Frequency')

title(str_title)

str_kp2_conf_int = [save_subfolder slash Param_name{1} '_' 'Confidence_Interval' '_' data_name '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

saveas(PL1,[str_kp2_conf_int,'.fig'], 'fig')

saveas(PL1,[str_kp2_conf_int,'.jpg'], 'jpg')

str_Best_Param_hist = [save_subfolder slash 'Best_Param_hist' '_' data_name '_' Param_name{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

xlswrite([str_Best_Param_hist, '.xlsx'],Ave_Prob_grid,'Sheet1','B1');

xlswrite([str_Best_Param_hist, '.xlsx'],kp2_bins,'Sheet1','A1');

col_header = {Param_name{1}};

row_header(1:4,1) = {'Mean', '1-sigma LB', '1-sigma UB', 'std_of_mean'};

str_FitResult = [save_subfolder slash 'Fit_outcome' '_' data_name '_' Param_name{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

xlswrite(str_FitResult,Fit_dat,'Sheet1','B2');

xlswrite(str_FitResult,col_header,'Sheet1','B1');

xlswrite(str_FitResult,row_header,'Sheet1','A2'); 


