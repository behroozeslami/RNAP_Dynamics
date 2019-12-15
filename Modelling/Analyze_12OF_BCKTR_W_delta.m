mainfoldername = 'Select_OF_12.5pN_1mM_NTP';

data_name = '12OF_BT';

grid_name = 'ML_grid_12OF_BT';

simul_name = 'Simul_DWTpdf_12OF_BT';

proc_name = 'Processivity_12OF_BT';

save_folder = 'Analysis_12OF_BT_W-delta';

subfolder_name = data_name;

t_strt = 0.05;

t_end = 10000;

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 4;

nboot = 100;

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

Param_name_1 = {'W'};

Param_name_2 = {'delta'};

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

str_simulgrid = [folder_grid slash simul_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

Simul_dat = load([str_simulgrid, '.mat'], '-mat');

Simul_grid = Simul_dat.Simul_grid;

str_Proc = [folder_grid slash proc_name '_' Param_name_1{1} '-' Param_name_2{1} '.mat'];

Proc_grid = load(str_Proc, '-mat');

Proc_grid = Proc_grid.Proc_grid;

str_DWT_dat_single = [output_foldername slash 'DWT_data' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

D = load(str_DWT_dat_single, '-mat');

DWT_Mat = D.DWT_dat; 

DWT_array = reshape(DWT_Mat,[1, numel(DWT_Mat)]);

DWT_array = DWT_array(DWT_array>0);

bootstat = bootstrp(nboot,@(DWT_array) DwellTimeHist_v3(DWT_array, t0, bins),DWT_array);

Exp_DWT_hist = mean(bootstat);
    
Err = std(bootstat);

binwidth = bins(2:length(bins)) - bins(1:length(bins)-1);

Simul_grid_P = Simul_grid(:, 3:length(bincenters)+2).*binwidth;

tail_P = sum(Simul_grid_P(:,bincenters>10),2);

[sorted_tail_P, idx] = sort(tail_P);

sorted_Simul_grid = Simul_grid(idx,:);

[sort_ML idx2] = sort(grid(:,3));

sort_grid = grid(idx2,:);

sort_grid = sort_grid(sort_grid(:,1)<11,:);

best_Wall = sort_grid(1,1);

best_delta = sort_grid(1,2);

best_ML = sort_grid(1,3);

best_ML_ourmodel = 4.0521e+03;

best_fit_data = [best_Wall, best_delta, best_ML, best_ML_ourmodel];

[sort_Proc, idx3] = sort(Proc_grid(:,3));

sort_Proc_grid = Proc_grid(idx3,:);

v = sorted_Simul_grid(:,1) == best_Wall & sorted_Simul_grid(:,2) == best_delta;

best_hist = sorted_Simul_grid(v, 3:length(bincenters)+2);

PL1 = figure;

ax = axes();

cmap = colormap(jet(length(sorted_Simul_grid)));

for n=1:length(Simul_grid)
    
    lw = 0.5;
    
    c = cmap(n,:);
    
    if ~(sorted_Simul_grid(n,1) == best_Wall && sorted_Simul_grid(n,2) == best_delta)
        
        dwt_hist = sorted_Simul_grid(n, 3:length(bincenters)+2);
    
        plot(bincenters, dwt_hist, '-', 'LineWidth', lw, 'Color',c)
    
        hold on
        
    end

end

plot(bincenters, best_hist, '-', 'LineWidth', 2, 'Color','k')

errorbar(bincenters,Exp_DWT_hist,Err,'kO','MarkerSize',6, 'MarkerFaceColor', 'k');

str_title = ['12OF' ' ' 'Ts = ' ' ' num2str((2*W+1)*t0) 's' ' ' 'Nw = ' num2str(Nw)];

title(str_title)

xlabel('Dwell Time (s)')

ylabel('PDF')

xlim([0.04, 500])

ylim([10^-5,2])

hold off

set(ax,'XScale','log','YScale','log')

PL2 = figure;

bar(Proc_grid(:,3))

xlabel('Parameters (Wall, \delta)')

ylabel('Processivity')

title('12OF')

max_Proc_grid = Proc_grid(Proc_grid(:,2) == 1,:);

min_Proc_grid = Proc_grid(Proc_grid(:,2) == 0,:);

str_Exp_Proc = [output_foldername slash 'Behrouz OF 12.5 pN values.xlsx'];

Exp_Proc_data = xlsread(str_Exp_Proc);

Exp_Proc = nanmean(Exp_Proc_data(:,1));

PL3 = figure;

semilogy(max_Proc_grid(:,1), max_Proc_grid(:,3), 'r')

hold on

semilogy(min_Proc_grid(:,1), min_Proc_grid(:,3), 'b')

line(xlim, [Exp_Proc, Exp_Proc], 'Color','k')

hold off

xlabel('Wall position')

ylabel('Processivity')

title('12OF')

legend('Max. Processivity (\delta = 1)', 'Min. Processivity (\delta = 0)', 'Exp. Processivity')

col_name_1 = {'Wall', 'delta'};

for j=1:length(bincenters)
    
    col_name_1{j+2} = ['bin_' num2str(j)];
    
end

str_simul_save = [save_subfolder slash 'Sorted_Simul_DWTpdf_12OF_BT' '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

xlswrite([str_simul_save,'.xlsx'],sorted_Simul_grid,'Sheet1','A2');

xlswrite([str_simul_save,'.xlsx'],col_name_1,'Sheet1','A1');

col_name_2 = {'bins', 'bin centers' 'DWT hist' 'Error'};

str_exp_save = [save_subfolder slash 'Exp_DWTpdf_12OF' '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

DWT_Exp_dat = [bins', [bincenters, 0]', [Exp_DWT_hist, 0]', [Err, 0]'];
    
DWT_Exp_dat(length(bins),2:4) = NaN([1,3]);

xlswrite([str_exp_save,'.xlsx'],DWT_Exp_dat,'Sheet1','A2');

xlswrite([str_exp_save,'.xlsx'],col_name_2,'Sheet1','A1');

str_exp_save = [save_subfolder slash 'Exp_DWTpdf_12OF' '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

col_name_3 = {'Wall', 'delta', 'ML score', 'ML score (our model)'};

str_fit = [save_subfolder slash 'BestFit_12OF_BT' '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

xlswrite([str_fit,'.xlsx'],best_fit_data,'Sheet1','A2');

xlswrite([str_fit,'.xlsx'],col_name_3,'Sheet1','A1');

col_name_4 = {'Wall', 'delta', 'ML score'};

str_ML = [save_subfolder slash 'Sorted_MLscore_12OF_BT' '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

xlswrite([str_ML,'.xlsx'],grid(idx2,:),'Sheet1','A2');

xlswrite([str_ML,'.xlsx'],col_name_4,'Sheet1','A1');

col_name_5 = {'Wall', 'delta', 'Processivity', 'Error'};

str_Proc_save = [save_subfolder slash 'Processivity_12OF_BT' '_' Param_name_1{1} '-' Param_name_2{1}];

xlswrite([str_Proc_save,'.xlsx'],Proc_grid,'Sheet1','A2');

xlswrite([str_Proc_save,'.xlsx'],col_name_5,'Sheet1','A1');

str_Proc_save = [save_subfolder slash 'Processivity_12OF_BT' '_' Param_name_1{1} '-' Param_name_2{1}];

xlswrite([str_Proc_save,'.xlsx'],Proc_grid,'Sheet1','A2');

xlswrite([str_Proc_save,'.xlsx'],col_name_5,'Sheet1','A1');

str_maxProc_save = [save_subfolder slash 'MaxProcessivity_12OF_BT' '_' Param_name_1{1} '-' Param_name_2{1}];

xlswrite([str_maxProc_save,'.xlsx'],max_Proc_grid,'Sheet1','A2');

xlswrite([str_maxProc_save,'.xlsx'],col_name_5,'Sheet1','A1');

str_minProc_save = [save_subfolder slash 'MinProcessivity_12OF_BT' '_' Param_name_1{1} '-' Param_name_2{1}];

xlswrite([str_minProc_save,'.xlsx'],min_Proc_grid,'Sheet1','A2');

xlswrite([str_minProc_save,'.xlsx'],col_name_5,'Sheet1','A1');

str_simul_save_2 = [save_subfolder slash 'Simul_DWTpdf_12OF_BT' '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

saveas(PL1,[str_simul_save_2,'.fig'], 'fig')

saveas(PL1,[str_simul_save_2,'.jpg'], 'jpg')

saveas(PL2,[str_Proc_save,'.fig'], 'fig')

saveas(PL2,[str_Proc_save,'.jpg'], 'jpg')

str_MaxMin_Proc_save = [save_subfolder slash 'Max-Min_Processivity_12OF_BT' '_' Param_name_1{1}];

saveas(PL3,[str_MaxMin_Proc_save,'.fig'], 'fig')

saveas(PL3,[str_MaxMin_Proc_save,'.jpg'], 'jpg')

