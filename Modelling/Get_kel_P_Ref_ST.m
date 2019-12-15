mainfoldername = 'Select_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

fit_index_1 = 1;

fit_index_2 = 2;

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

save_folder = 'Analysis_1mM_SingleTraces_kel-P';

subfolder_name = '1mM_indv';

t_strt = 0.05;

t_end = 10000;

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Nw = 10;

[Param_dat, txt] = xlsread('Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Param_name_1 = txt(1,fit_index_1+1);

Param_name_2 = txt(1,fit_index_2+1);

save_subfolder = [output_foldername '/' save_folder '/' subfolder_name];

trace_indices = load([output_foldername slash 'Trace_indices.mat']);

trace_indices = trace_indices.selected_trace_number;

data_kel_P = zeros([length(trace_indices), 2]);

for j=1:length(trace_indices)
    
    data_name = ['1mM_' num2str(trace_indices(j))];
    
    str_FitResult = [save_subfolder slash 'Fit_outcome' '_' data_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];
    
    data = xlsread(str_FitResult);
    
    data_kel_P(j,:) = data(1,:);
    
end

str_out = [output_foldername '/' 'Fit_outcome_1mM_kel_P.txt'];

save(str_out, 'data_kel_P', '-ascii')


