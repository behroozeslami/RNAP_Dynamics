W = 12;

k = 0;

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

file_str_txt = [output_foldername slash 'Trace_data' '.txt'];

dat = load(file_str_txt, '-ascii');

l0_WLC = dat(1,4);

str_trace_index = [DWT_folder slash 'Trace_indices' '.mat'];

str_V = [DWT_folder slash 'EndtoEnd_V' '.mat'];

trace_indices = load(str_trace_index);

trace_indices = trace_indices.selected_trace_number;

L = length(trace_indices);

V_array = zeros([1,L]);

for j=1:L
    
    filename = ['RNAP_' num2str(trace_indices(j)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = trace_dat(1,10);
    
    factor = l0/l0_WLC;
    
    %factor = 1;

    x_array = factor*trace_dat(:,8);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);
    
    Max_I = length(xs_array);
    
    delta_x = xs_array(Max_I)-xs_array(1);
    
    delta_T = t_array(Max_I)-t_array(1);
    
    V = delta_x/delta_T;
    
    V_array(j) = V;
    
end

save(str_V, 'V_array', '-mat')