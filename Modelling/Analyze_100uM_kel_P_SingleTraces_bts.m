function Analyze_100uM_kel_P_SingleTraces_bts(data_name, subfolder_name, selected_trace_indices) 

mainfoldername = 'NEW_AF_7.5pN_100uM_NTP';

grid_name = 'ML_grid_100uM';

simul_name = 'Simul_DWTpdf_100uM';

save_folder = 'Analysis_100uM_SingleTraces_kel-P';

fit_index_1 = 1;

fit_index_2 = 2;

param_num = 2;

t_strt = 0.05;

t_end = 10000;

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 4;

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

trace_indices = load([output_foldername slash 'Trace_indices_withBS.mat']);

trace_indices = trace_indices.selected_trace_number;

Unique_trace_indices = load([output_foldername slash 'Trace_indices.mat']);

Unique_trace_indices = Unique_trace_indices.selected_trace_number;

subset_cols =ismember(trace_indices, selected_trace_indices);

ndat = length(trace_indices);

[Param_dat, txt] = xlsread('Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Params = Param_dat(1,:);

Param_name_1 = txt(1,fit_index_1+1);

Param_name_2 = txt(1,fit_index_2+1);

folder_grid = [output_foldername slash save_folder];

str_grid = [folder_grid slash grid_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

grid_dat = load(str_grid,'-mat');

grid = grid_dat.grid;

str_bincenters = [output_foldername '/' 'DWT_pdf_bincenters' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];
 
bincenters = load(str_bincenters, '-mat');

bincenters = bincenters.centers;

str_bins = [output_foldername '/' 'DWT_pdf_bins' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];
 
bins = load(str_bins, '-mat');

bins = bins.bins;

str_simulgrid = [folder_grid slash simul_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];

Simul_dat = load(str_simulgrid, '-mat');

Simul_grid = Simul_dat.Simul_grid;

S = size(Simul_grid);

nbins = S(2)-2;

str_datlength = [output_foldername '/' 'DWT_data_length_withBS' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

DL = load(str_datlength, '-mat');

dat_length = DL.dat_length;

str_unique_datlength = [output_foldername '/' 'DWT_data_length' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

DL = load(str_unique_datlength, '-mat');

str_DWT_dat_single = [output_foldername slash 'DWT_data' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

D = load(str_DWT_dat_single, '-mat');

DWT_Mat = D.DWT_dat; 

Unique_dat_length = DL.dat_length;

temp_grid = grid(:,param_num+1:ndat+param_num);

temp_grid = temp_grid(:,subset_cols);

grid = [grid(:,1:param_num), temp_grid];

dat_length = dat_length(subset_cols);

trace_indices = trace_indices(subset_cols);

Unique_dat_length = Unique_dat_length(ismember(Unique_trace_indices, selected_trace_indices));

DWT_Mat = DWT_Mat(ismember(Unique_trace_indices, selected_trace_indices),:);

DWT_array = reshape(DWT_Mat,[1, numel(DWT_Mat)]);

DWT_array = DWT_array(DWT_array>0);

Unique_trace_indices = intersect(Unique_trace_indices,selected_trace_indices); 

ndat = length(dat_length);

[Min_LL, best_indices] = min(grid(:,param_num+1:ndat+param_num));

gSize = size(grid);

points_num = gSize(1);

bestparam_grid = zeros([points_num, ndat]);

for n=1:ndat
    
    ind = best_indices(n);
    
    bestparam_grid(ind,n) = 1;
    
end

bestparam_agg_grid = cell2mat(arrayfun(@(I) sum((bestparam_grid*transpose(ismember(trace_indices, I))), 2), Unique_trace_indices, 'UniformOutput', false));

Prob_grid = cell2mat(arrayfun(@(j) remove_outliers_and_average(bestparam_agg_grid(:,j), param_num, grid), 1:length(Unique_trace_indices), 'UniformOutput', false));

dat_length_normal = Unique_dat_length/sum(Unique_dat_length);

Ave_Prob_grid = Prob_grid*transpose(dat_length_normal);

kel_bins = unique(grid(:,1));

P_bins = unique(grid(:,2));

delta_kel = kel_bins(2)-kel_bins(1);

delta_P = P_bins(2)-P_bins(1);

P_Prob_grid = cell2mat(arrayfun(@(P) sum(Ave_Prob_grid(grid(:,2)==P,1)), P_bins,'UniformOutput',0));

kel_Prob_grid = cell2mat(arrayfun(@(kel) sum(Ave_Prob_grid(grid(:,1)==kel,1)), kel_bins,'UniformOutput',0));

[Best_P, Mean_P, LB_P, UB_P, std_P] = Conf_int(P_bins,P_Prob_grid,alpha);

[Best_kel, Mean_kel, LB_kel, UB_kel, std_kel] = Conf_int(kel_bins,kel_Prob_grid,alpha);

Fit_dat = [[Mean_kel, LB_kel, UB_kel, std_kel/sqrt(length(selected_trace_indices))]', [Mean_P, LB_P, UB_P, std_P/sqrt(length(selected_trace_indices))]'];

str_title = strrep(data_name, '_', ' - ');

PL1 = figure(1);

bar(kel_bins,kel_Prob_grid,'r')

hold on

line([LB_kel-0.5*delta_kel,LB_kel-0.5*delta_kel],ylim,'LineWidth', 2, 'Color', 'b')

line([UB_kel+0.5*delta_kel,UB_kel+0.5*delta_kel],ylim,'LineWidth', 2, 'Color', 'b')

line([Mean_kel,Mean_kel],ylim,'LineWidth', 2, 'Color', 'g')

hold off

xlabel('Optimal kel (bp/s)')

ylabel('Frequency')

title(str_title)

PL2 = figure(2);

bar(P_bins,P_Prob_grid,'r')

hold on

line([LB_P-0.5*delta_P,LB_P-0.5*delta_P],ylim,'LineWidth', 2, 'Color', 'b')

line([UB_P+0.5*delta_P,UB_P+0.5*delta_P],ylim,'LineWidth', 2, 'Color', 'b')

line([Mean_P,Mean_P],ylim,'LineWidth', 2, 'Color', 'g')

hold off

xlabel('Optimal P')

ylabel('Frequency')

title(str_title)

Ave_Prob_mat = zeros([length(P_bins), length(kel_bins)]);

for n=1:length(P_bins)
    
    for m=1:length(kel_bins)
        
        P = P_bins(n);
        
        kel = kel_bins(m);
        
        Ave_Prob_mat(n,m) = Ave_Prob_grid(((grid(:,2)==P)&(grid(:,1)==kel)),1);
        
    end
    
end

PL3 = figure(3);

surf(kel_bins, P_bins, Ave_Prob_mat);

view([0,0,1])

colormap('jet')

shading interp

colorbar

xlabel('kb (bp/s)')

ylabel('P')

zlabel('Probability')

title(str_title)

Simul_DWT_hist_combined = transpose(Ave_Prob_grid)*Simul_grid(:, param_num+1:nbins+param_num);

[Min,Min_Ind] = min((Mean_kel-grid(:,1)).^2 + (Mean_P-grid(:,2)).^2);

Simul_DWT_hist_mean = Simul_grid(Min_Ind, param_num+1:nbins+param_num);

nboot = 100;

bootstat = bootstrp(nboot,@(DWT_array) DwellTimeHist_v3(DWT_array, t0, bins),DWT_array);

Exp_DWT_hist = mean(bootstat);
    
Err = std(bootstat);

PL4 = figure(4);

ax = axes();

errorbar(bincenters,Exp_DWT_hist,Err,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

hold on

loglog(bincenters, Simul_DWT_hist_combined, 'k-', 'LineWidth', 2)

loglog(bincenters, Simul_DWT_hist_mean, 'g--', 'LineWidth', 2)

hold off

legend('Data' , 'Best Fit (combined)', 'Best Fit (mean)')

str_title = [str_title ' ' 'Ts = ' ' ' num2str((2*W+1)*t0) 's' ' ' 'Nw = ' num2str(Nw)];

title(str_title)

xlabel('Dwell Time (s)')

ylabel('PDF')

set(ax,'XScale','log','YScale','log')

str_kel_conf_int = [save_subfolder slash Param_name_2{1} '_' 'Confidence_Interval' '_' data_name '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

str_P_conf_int = [save_subfolder slash Param_name_1{1} '_' 'Confidence_Interval' '_' data_name '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

saveas(PL1,[str_kel_conf_int,'.fig'], 'fig')

saveas(PL1,[str_kel_conf_int,'.jpg'], 'jpg')

saveas(PL2,[str_P_conf_int,'.fig'], 'fig')

saveas(PL2,[str_P_conf_int,'.jpg'], 'jpg')

str_Best_Param_hist = [save_subfolder slash 'Best_Param_2Dhist' '_' data_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

saveas(PL3,[str_Best_Param_hist,'.fig'], 'fig')

saveas(PL3,[str_Best_Param_hist,'.jpg'], 'jpg')

xlswrite([str_Best_Param_hist, '.xlsx'],Ave_Prob_mat,'Sheet1','C3');

xlswrite([str_Best_Param_hist, '.xlsx'],P_bins,'Sheet1','A3');

xlswrite([str_Best_Param_hist, '.xlsx'],kel_bins','Sheet1','C1');

xlswrite([str_Best_Param_hist, '.xlsx'],Param_name_1','Sheet1','B1');

xlswrite([str_Best_Param_hist, '.xlsx'],Param_name_2,'Sheet1','A2');

col_header = {Param_name_1{1}, Param_name_2{1}};

row_header(1:4,1) = {'Mean', '1-sigma LB', '1-sigma UB', 'std_of_mean'};

str_FitResult = [save_subfolder slash 'Fit_outcome' '_' data_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

xlswrite(str_FitResult,Fit_dat,'Sheet1','B2');

xlswrite(str_FitResult,col_header,'Sheet1','B1');

xlswrite(str_FitResult,row_header,'Sheet1','A2'); 

str_DWT_Simul = [save_subfolder slash 'Best_fit_DWTpdf' '_' data_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

saveas(PL4,[str_DWT_Simul,'.fig'], 'fig')

saveas(PL4,[str_DWT_Simul,'.jpg'], 'jpg')

DWT_Simul_dat = [bins', [bincenters, 0]', [Exp_DWT_hist, 0]', [Err, 0]', [Simul_DWT_hist_combined, 0]', [Simul_DWT_hist_mean, 0]'];
    
DWT_Simul_dat(length(bins),2:6) = NaN([1,5]);

header = {'Bin edges', 'Bin centers', 'Exp DWT pdf', 'Errors', 'Best fit (combined)', 'Best fit (mean)'};
    
xlswrite([str_DWT_Simul,'.xlsx'],DWT_Simul_dat,'Sheet1','A2');

xlswrite([str_DWT_Simul,'.xlsx'],header,'Sheet1','A1');

close all



