function Analyze_1uM_kel_P_Combined_bts(data_name, subfolder_name, selected_trace_indices) 

mainfoldername = 'NEW_AF_7.5pN_1uM_NTP';

grid_name = 'ML_grid_1uM';

simul_name = 'Simul_DWTpdf_1uM';

save_folder = 'Analysis_1uM_SingleTraces_kel-P';

fit_index_1 = 1;

fit_index_2 = 2;

param_num = 2;

nboot = 50;

t_strt = 0.05;

t_end = 10000;

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 4;

alpha = 0.136; % 1-sigma confidence interval

girid_kel_ub = 50;

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

[Param_dat, txt] = xlsread('Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Params = Param_dat(1,:);

Param_name_1 = txt(1,fit_index_1+1);

Param_name_2 = txt(1,fit_index_2+1);

folder_grid = [output_foldername slash save_folder];

str_grid = [folder_grid slash grid_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

grid_dat = load(str_grid,'-mat');

grid = grid_dat.grid;

sum_LL_grid = Calculate_LogLikelihood_Combine_from_SingleTrace(grid,selected_trace_indices,trace_indices,param_num,nboot);

kel_bins = unique(grid(:,1));

P_bins = unique(grid(:,2));

sum_LL_kel = zeros([length(kel_bins),nboot]);

for i=1:length(kel_bins)
    
    sum_LL_kel(i,:) = min(sum_LL_grid(grid(:,1)==kel_bins(i),:));
    
end

delta_sum_LL_kel = cell2mat(transpose(arrayfun(@(row) sum_LL_kel(row,:)-min(sum_LL_kel),1:length(kel_bins), 'UniformOutput', 0)));

Mean_LL_kel = mean(delta_sum_LL_kel,2);

Mean_LL_kel = Mean_LL_kel - min(Mean_LL_kel);

sigma_LL_kel = std(delta_sum_LL_kel,0,2);

inds = 1:length(kel_bins);

inds = inds(Mean_LL_kel<11);

ind = min(inds);

kel_LB = kel_bins(ind);

%{
func = @(P2,P4,x0,x) P2*(x-x0).^2+P4*(x-x0).^4;

fit_LL = fittype(func);

x = kel_bins;

y = Mean_LL_kel;

Smooth_LL = fit(x,y,fit_LL, 'Start', [1,1,1], 'Lower',[0,0,0]);
%}

str_title = strrep(data_name, '_', ' - ');

figure(1)

%plot(Smooth_LL,'b-');

hold on

PL1 = errorbar(kel_bins, Mean_LL_kel,sigma_LL_kel,'O-', 'Color','r');

plot([kel_LB,kel_LB],ylim)

hold off

xlabel(Param_name_1)

ylabel('minus LogLikelihood')

title(str_title)

sum_LL_P = zeros([length(P_bins),nboot]);

for i=1:length(P_bins)
    
    sum_LL_P(i,:) = min(sum_LL_grid(grid(:,2)==P_bins(i),:));
    
end

delta_sum_LL_P = cell2mat(transpose(arrayfun(@(row) sum_LL_P(row,:)-min(sum_LL_P),1:length(P_bins), 'UniformOutput', 0)));

Mean_LL_P = mean(delta_sum_LL_P,2);

sigma_LL_P = std(delta_sum_LL_P,0,2);

figure(2)

PL2 = errorbar(P_bins, Mean_LL_P,sigma_LL_P,'O-', 'Color','r');

xlabel(Param_name_2)

ylabel('minus LogLikelihood')

title(str_title)

str_kel_LL = [save_subfolder slash Param_name_1{1} '_' 'minus_LogLikelihood' '_' data_name '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

str_P_LL = [save_subfolder slash Param_name_2{1} '_' 'minus_LogLikelihood' '_' data_name '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

saveas(PL1,[str_kel_LL,'.fig'], 'fig')

saveas(PL1,[str_kel_LL,'.jpg'], 'jpg')

saveas(PL2,[str_P_LL,'.fig'], 'fig')

saveas(PL2,[str_P_LL,'.jpg'], 'jpg')

bound_LL_grid = sum_LL_grid(grid(:,1)<=girid_kel_ub,:);

[Min_LL, best_indices] = min(bound_LL_grid);

best_kel = grid(best_indices,1);

best_P = grid(best_indices,2);

[Mean_P, LB_P, UB_P, std_P] = Conf_int_4_combined(best_P,alpha);

[Mean_kel, LB_kel, UB_kel, std_kel] = Conf_int_4_combined(best_kel,alpha);

Fit_dat = [[Mean_kel, LB_kel, UB_kel, std_kel, kel_LB]', [Mean_P, LB_P, UB_P, std_P, NaN(1)]'];

col_header = {Param_name_1{1}, Param_name_2{1}};

row_header(1:5,1) = {'Mean', '1-sigma LB', '1-sigma UB', 'std', 'kel_LB_LL10'};

str_FitResult = [save_subfolder slash 'Fit_outcome' '_' data_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '_', 'kel_UB = ' num2str(girid_kel_ub) '.xlsx'];

xlswrite(str_FitResult,Fit_dat,'Sheet1','B2');

xlswrite(str_FitResult,col_header,'Sheet1','B1');

xlswrite(str_FitResult,row_header,'Sheet1','A2'); 

[Min,Min_Ind] = min((Mean_kel-grid(:,1)).^2 + (Mean_P-grid(:,2)).^2);

str_simulgrid = [folder_grid slash simul_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];

Simul_dat = load(str_simulgrid, '-mat');

Simul_grid = Simul_dat.Simul_grid;

S = size(Simul_grid);

nbins = S(2)-2;

Simul_DWT_hist_mean = Simul_grid(Min_Ind, param_num+1:nbins+param_num);

str_DWT_dat_single = [output_foldername slash 'DWT_data' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

D = load(str_DWT_dat_single, '-mat');

DWT_Mat = D.DWT_dat; 

Unique_trace_indices = load([output_foldername slash 'Trace_indices.mat']);

Unique_trace_indices = Unique_trace_indices.selected_trace_number;

DWT_Mat = DWT_Mat(ismember(Unique_trace_indices, selected_trace_indices),:);

DWT_array = reshape(DWT_Mat,[1, numel(DWT_Mat)]);

DWT_array = DWT_array(DWT_array>0);

str_bincenters = [output_foldername '/' 'DWT_pdf_bincenters' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];
 
bincenters = load(str_bincenters, '-mat');

bincenters = bincenters.centers;

str_bins = [output_foldername '/' 'DWT_pdf_bins' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];
 
bins = load(str_bins, '-mat');

bins = bins.bins;

bootstat = bootstrp(nboot,@(DWT_array) DwellTimeHist_v3(DWT_array, t0, bins),DWT_array);

Exp_DWT_hist = mean(bootstat);
    
Err = std(bootstat);

PL4 = figure(4);

ax = axes();

errorbar(bincenters,Exp_DWT_hist,Err,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

hold on

loglog(bincenters, Simul_DWT_hist_mean, 'g-', 'LineWidth', 2)

hold off

legend('Data' , 'Best Fit')

str_title = [str_title ' ' 'Ts = ' ' ' num2str((2*W+1)*t0) 's' ' ' 'Nw = ' num2str(Nw)];

title(str_title)

xlabel('Dwell Time (s)')

ylabel('PDF')

set(ax,'XScale','log','YScale','log')

str_DWT_Simul = [save_subfolder slash 'Best_fit_DWTpdf' '_' data_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

saveas(PL4,[str_DWT_Simul,'.fig'], 'fig')

saveas(PL4,[str_DWT_Simul,'.jpg'], 'jpg')

DWT_Simul_dat = [bins', [bincenters, 0]', [Exp_DWT_hist, 0]', [Err, 0]', [Simul_DWT_hist_mean, 0]'];
    
DWT_Simul_dat(length(bins),2:5) = NaN([1,4]);

header = {'Bin edges', 'Bin centers', 'Exp DWT pdf', 'Errors', 'Best fit'};
    
xlswrite([str_DWT_Simul,'.xlsx'],DWT_Simul_dat,'Sheet1','A2');

xlswrite([str_DWT_Simul,'.xlsx'],header,'Sheet1','A1');

close all
