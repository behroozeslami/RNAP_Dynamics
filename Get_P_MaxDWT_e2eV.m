W = 12;

k = 0;

Nw = 4;

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0 = 1/cam_freq;

t_strt = 0.05;

t_end = 10000;

T_cut = (2*W+1)*t0;

mainfoldername = 'NEW_AF_7.5pN_1mM_ITP_100uM_NTP';

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

str_trace_index = [DWT_folder slash 'Trace_indices' '.mat'];

trace_indices = load(str_trace_index);

trace_indices = trace_indices.selected_trace_number;

str_DWT_length = [DWT_folder slash 'DWT_data_length' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

dat_length = load(str_DWT_length);

DWT_length = dat_length.dat_length;

str_DWT_dat = [DWT_folder slash 'DWT_data' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

DWTs = load(str_DWT_dat);

DWTs = DWTs.DWT_dat;

str_V = [DWT_folder slash 'EndtoEnd_V' '.mat'];

V = load(str_V);

V = V.V_array;

P = arrayfun(@(row) sum(DWTs(row, :)>T_cut)/DWT_length(row), 1:length(trace_indices));

p_bins = 0.01:0.01:ceil(100*max(P))/100;

P_hist = hist(P, p_bins);

P_hist = P_hist/sum(P_hist);

figure(1)

bar(p_bins, P_hist)

Max_DWT = arrayfun(@(row) max(DWTs(row, :)), 1:length(trace_indices));

Max_DWT_bins = 0:10:ceil(10*max(Max_DWT))/10;

MAX_DWT_hist = hist(Max_DWT, Max_DWT_bins);

MAX_DWT_hist = MAX_DWT_hist/sum(MAX_DWT_hist);

figure(2)

bar(Max_DWT_bins, MAX_DWT_hist)

figure(3)

semilogy(P,Max_DWT, 'r.')

V_bins = 0:0.1:ceil(max(V)*10)/10;

V_hist = hist(V,V_bins);

V_hist = V_hist/sum(V_hist);

figure(4)

bar(V_bins, V_hist)

figure(5)

semilogy(V,Max_DWT, 'r.')

figure(6)

plot(V, P, 'r.')

