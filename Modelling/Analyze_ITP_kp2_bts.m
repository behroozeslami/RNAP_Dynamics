function Analyze_ITP_kp2_bts(data_name) 

mainfoldername = 'NEW_AF_7.5pN_1mM_ITP_100uM_NTP';

grid_name = 'ML_grid_ITP';

simul_name = 'Simul_DWTpdf_ITP';

save_folder = 'Analysis_ITP_kp2';

fit_index_1 = 6;

param_num = 1;

t_strt = 0.05;

t_end = 30;

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

subfolder_name = data_name;

save_subfolder = [output_foldername '/' save_folder '/' subfolder_name];

if ~(exist(save_subfolder,'dir')==7)
    
    mkdir(save_subfolder)
    
end

[Param_dat, txt] = xlsread('Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Params = Param_dat(1,:);

Param_name_1 = txt(1,fit_index_1+1);

folder_grid = [output_foldername slash save_folder];

str_grid = [folder_grid slash grid_name '_' Param_name_1{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

grid_dat = load(str_grid,'-mat');

grid = grid_dat.grid;

str_bincenters = [output_foldername '/' 'DWT_pdf_bincenters' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];
 
bincenters = load(str_bincenters, '-mat');

bincenters = bincenters.centers;

str_bins = [output_foldername '/' 'DWT_pdf_bins' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];
 
bins = load(str_bins, '-mat');

bins = bins.bins;

str_simulgrid = [folder_grid slash simul_name '_' Param_name_1{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];

Simul_dat = load(str_simulgrid, '-mat');

Simul_grid = Simul_dat.Simul_grid;

S = size(Simul_grid);

nbins = S(2)-param_num;

str_datlength = [output_foldername '/' 'DWT_data_length_withBS' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

DL = load(str_datlength, '-mat');

dat_length = DL.dat_length;

str_DWT_dat = [output_foldername slash 'DWT_data' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(10000) '.mat'];

D = load(str_DWT_dat, '-mat');

DWT_Mat = D.DWT_dat; 

DWT_array = reshape(DWT_Mat,[1, numel(DWT_Mat)]);

DWT_array = DWT_array(DWT_array>0);

ndat = length(dat_length);

[Min_LL, best_indices] = min(grid(:,param_num+1:ndat+param_num));

gSize = size(grid);

points_num = gSize(1);

bestparam_grid = zeros([points_num, ndat]);

for n=1:ndat
    
    ind = best_indices(n);
    
    bestparam_grid(ind,n) = 1;
    
end

Ave_Prob_grid = bestparam_grid*ones([ndat,1])/ndat;

kp2_bins = unique(grid(:,1));

delta_kp2 = kp2_bins(2)-kp2_bins(1);

kp2_Prob_grid = cell2mat(arrayfun(@(kp2) sum(Ave_Prob_grid(grid(:,1)==kp2,1)), kp2_bins,'UniformOutput',0));

[Best_kp2, Mean_kp2, LB_kp2, UB_kp2, std_kp2] = Conf_int(kp2_bins,kp2_Prob_grid,alpha);

if abs(LB_kp2 - Best_kp2) < delta_kp2
    
    LB_kp2 = LB_kp2 - delta_kp2;
    
end

if abs(UB_kp2 - Best_kp2) < delta_kp2
    
    UB_kp2 = UB_kp2 + delta_kp2;
    
end

if std_kp2 < delta_kp2
    
    std_kp2 = delta_kp2;
    
end
    

Fit_dat = [Mean_kp2, LB_kp2, UB_kp2, std_kp2]';

str_title = strrep(data_name, '_', ' - ');

PL1 = figure(1);

bar(kp2_bins,kp2_Prob_grid,'r')

hold on

line([LB_kp2-0.5*delta_kp2,LB_kp2-0.5*delta_kp2],ylim,'LineWidth', 2, 'Color', 'b')

line([UB_kp2+0.5*delta_kp2,UB_kp2+0.5*delta_kp2],ylim,'LineWidth', 2, 'Color', 'b')

line([Mean_kp2,Mean_kp2],ylim,'LineWidth', 2, 'Color', 'g')

hold off

xlabel('Optimal kp2 (1/s)')

ylabel('Frequency')

title(str_title)

Simul_DWT_hist_combined = transpose(Ave_Prob_grid)*Simul_grid(:, param_num+1:nbins+param_num);

[Min,Min_Ind] = min((Mean_kp2-grid(:,1)).^2);

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

str_kp2_conf_int = [save_subfolder slash Param_name_1{1} '_' 'Confidence_Interval' '_' data_name '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

saveas(PL1,[str_kp2_conf_int,'.fig'], 'fig')

saveas(PL1,[str_kp2_conf_int,'.jpg'], 'jpg')

col_header = {Param_name_1{1}};

row_header(1:4,1) = {'Mean', '1-sigma LB', '1-sigma UB', 'std_of_mean'};

str_FitResult = [save_subfolder slash 'Fit_outcome' '_' data_name '_' Param_name_1{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

xlswrite(str_FitResult,Fit_dat,'Sheet1','B2');

xlswrite(str_FitResult,col_header,'Sheet1','B1');

xlswrite(str_FitResult,row_header,'Sheet1','A2'); 

str_DWT_Simul = [save_subfolder slash 'Best_fit_DWTpdf' '_' data_name '_' Param_name_1{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

saveas(PL4,[str_DWT_Simul,'.fig'], 'fig')

saveas(PL4,[str_DWT_Simul,'.jpg'], 'jpg')

DWT_Simul_dat = [bins', [bincenters, 0]', [Exp_DWT_hist, 0]', [Err, 0]', [Simul_DWT_hist_combined, 0]', [Simul_DWT_hist_mean, 0]'];
    
DWT_Simul_dat(length(bins),2:6) = NaN([1,5]);

header = {'Bin edges', 'Bin centers', 'Exp DWT pdf', 'Errors', 'Best fit (combined)', 'Best fit (mean)'};
    
xlswrite([str_DWT_Simul,'.xlsx'],DWT_Simul_dat,'Sheet1','A2');

xlswrite([str_DWT_Simul,'.xlsx'],header,'Sheet1','A1');

close all



