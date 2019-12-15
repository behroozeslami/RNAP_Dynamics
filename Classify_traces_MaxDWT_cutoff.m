W = 24;

k = 0;

Nw = 4;

cam_freq = 25;
%%%%%% IMPORTANT
cam_freq = 50;

t0 = 1/cam_freq;

t_strt = 0.05;

t_end = 10000;

alpha = 0.1;

mainfoldername = 'OF4_9pN_1mM_NTP_GreA';

data_name = 'GreA';

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

select_folder = [DWT_folder slash 'Selection_maxDWT_cutoff=' num2str(alpha)];

if ~(exist(select_folder,'dir')==7)
    
    mkdir(select_folder)
    
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

Max_DWT = arrayfun(@(row) max(DWTs(row, :)), 1:length(trace_indices));

cutoff = quantile(Max_DWT,1-alpha);

selected_trace_indices = trace_indices(Max_DWT<cutoff);

selected_Max_DWT = Max_DWT(Max_DWT<cutoff);

Max_DWT_bins = 0:10:ceil(10*max(Max_DWT))/10;

colors = {'r','b'};

PL = histogram(Max_DWT, Max_DWT_bins, 'EdgeColor','k', 'FaceColor', 'b');

hold on

histogram(selected_Max_DWT, Max_DWT_bins, 'EdgeColor','k', 'FaceColor', 'r');

hold off

xlabel('Max. Dwell Time (s)')

ylabel('count')

title(data_name)

str_fig = [select_folder slash data_name, '_' 'MaxDWT_cutoff=' num2str(alpha)];

saveas(PL, [str_fig, '.bmp'], 'bmp')

saveas(PL, [str_fig, '.fig'], 'fig')

class_name = [select_folder slash data_name, '_', 'cotoff=', num2str(alpha) '.txt'];

save(class_name, 'selected_trace_indices', '-ascii')