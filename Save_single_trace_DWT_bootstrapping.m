W = 12;

k = 0;

Nw = 10;

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0 = 1/cam_freq;

t_strt = 0.05;

t_end = 10000;

nboot = 50;

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

DWT_folder = ['Modelling' slash output_foldername];

if ~(exist(DWT_folder,'dir')==7)
    
    mkdir(DWT_folder)
    
end

directory = ['..' slash 'Richard_Data' slash mainfoldername slash subfolder];

if strcmp('',subfolder)
    
    directory = ['..' slash 'Richard_Data' slash mainfoldername];
    
end

file_str_txt = [output_foldername slash 'Trace_data' '.txt'];

dat = load(file_str_txt, '-ascii');

trace_number = transpose(dat(:,1));

seidel_array = transpose(dat(:,2));

l0_array = transpose(dat(:,3));

l0_WLC = dat(1,4);

select_vec = (l0_array >= l0_min) & (l0_array <= l0_max) & (seidel_array >= Seidel_min) & (seidel_array <= Seidel_max);

selected_trace_number = floor(trace_number(select_vec));

selected_l0_array = l0_array(select_vec);

High_V_Blacklist = load([output_foldername '/' 'Blacklist_High_Velocity.txt'], '-ascii');

Low_V_Blacklist = load([output_foldername '/' 'Blacklist_Low_Velocity.txt'], '-ascii');

Blacklist = load([output_foldername '/' 'Blacklist.txt'], '-ascii');

Blacklist = [Blacklist; High_V_Blacklist; Low_V_Blacklist];

index_2_keep = Exclude_Blacklist(selected_trace_number, Blacklist);

selected_trace_number = selected_trace_number(index_2_keep);

selected_l0_array = selected_l0_array(index_2_keep);

str_DWT_length = [DWT_folder slash 'DWT_data_length_withBS' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

str_DWT_dat = [DWT_folder slash 'DWT_data_withBS' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

str_trace_index = [DWT_folder slash 'Trace_indices_withBS' '.mat'];

DWT_folder_0 = [output_foldername slash 'selected_DWT_pdf_Combined' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

str_DWT_0 = [DWT_folder_0 slash 'DWT_pdf_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

DWT_hist_dat_0 = load([str_DWT_0, '.txt'], '-ascii');

bins = transpose(DWT_hist_dat_0(:,1));

centers = transpose(DWT_hist_dat_0(1:length(bins)-1,2));

L = length(selected_trace_number);

dat_length = zeros([1,L]);

for i=1:L
    
    filename = ['RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = selected_l0_array(i);
    
    factor = l0/l0_WLC;
    
    %factor = 1;

    x_array = factor*trace_dat(:,8);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);

    [Tc_array]=find_crossing_time(t_array, xs_array,Nw);

    [DWT_array] = find_Dwell_time(Tc_array);
    
    DWT_array = DWT_array(DWT_array > t_strt & DWT_array < t_end);
        
    dat_length(i) = length(DWT_array);

end

DWT_dat = zeros([L, max(dat_length)]);

DWT_hist_dat = zeros([L, length(bins)-1]);

DWT_hist_err = zeros([L, length(bins)-1]);

for i=1:L
    
    filename = ['RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = selected_l0_array(i);
    
    factor = l0/l0_WLC;
    
    %factor = 1;

    x_array = factor*trace_dat(:,8);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);

    [Tc_array]=find_crossing_time(t_array, xs_array,Nw);

    [DWT_array] = find_Dwell_time(Tc_array);
    
    DWT_array = DWT_array(DWT_array > t_strt & DWT_array < t_end);
        
    DWT_dat(i,1:dat_length(i)) = DWT_array;

end

DWT_dat_BS = cell2mat(arrayfun(@(row) bootstrap_DWT(row, DWT_dat, dat_length, nboot), 1:L, 'UniformOutput', false)');

dat_length_BS = cell2mat(arrayfun(@(x) x*ones([1,nboot]), dat_length, 'UniformOutput', false));

selected_trace_number_BS = cell2mat(arrayfun(@(x) x*ones([1,nboot]), selected_trace_number, 'UniformOutput', false));

DWT_dat = DWT_dat_BS;

dat_length = dat_length_BS;

selected_trace_number = selected_trace_number_BS;

save(str_DWT_dat, 'DWT_dat', '-mat')

save(str_DWT_length, 'dat_length', '-mat')

save(str_trace_index, 'selected_trace_number', '-mat');
