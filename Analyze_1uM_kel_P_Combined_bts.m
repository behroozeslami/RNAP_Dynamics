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

[Min_LL, best_indices] = min(sum_LL_grid);

best_kel = grid(best_indices,1);

best_P = grid(best_indices,2);

kel_bins = unique(grid(:,1));

P_bins = unique(grid(:,2));

sum_LL_kel = zeros([length(kel_bins),nboot]);

for i=1:length(kel_bins)
    
    sum_LL_kel(i,:) = min(sum_LL_grid(grid(:,1)==kel_bins(i),:));
    
end