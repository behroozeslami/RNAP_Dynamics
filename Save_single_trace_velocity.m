W = 12;

k = 0;

Nw = 10;

Pause_threshold_V = 100.64;

T = (2*W+1);

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0 = 1/cam_freq;

Vmax_fit = 30;

Vmin_fit = -20;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'Export_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

V_folder = ['Modelling' slash output_foldername];

if ~(exist(V_folder,'dir')==7)
    
    mkdir(V_folder)
    
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

str_V_length = [V_folder slash 'V_data_length' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T)  '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '_' 'limits=' num2str(Vmin_fit) '-' num2str(Vmax_fit) '.mat'];

str_V_dat = [V_folder slash 'V_data' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '_' 'limits=' num2str(Vmin_fit) '-' num2str(Vmax_fit) '.mat'];

L = length(selected_trace_number);

dat_length = zeros([1,L]);

for i=1:L
       
    filename = ['RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = selected_l0_array(i);
    
    factor = l0/l0_WLC;

    x_array = factor*trace_dat(:,8);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);
    
    len = length(xs_array);
    
    Dist_array = transpose(xs_array(T+1:len) - xs_array(1:len-T));
    
    [Tc_array] = find_crossing_time(t_array, xs_array,Nw);
    
    [DWT_array]=find_Dwell_time(Tc_array);
    
    DWT_indices = 1:length(DWT_array);
    
    pause_start_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V));
    
    pause_end_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V)+1);
    
    pause_start_indices = round(pause_start_times/t0 +1);
    
    pause_end_indices = round(pause_end_times/t0 +1);
    
    pause_start_indices = pause_start_indices-(T+1);
    
    pause_start_indices(pause_start_indices < 1) = 1;
    
    bool = true([1, len]);
    
    for j=1:length(pause_start_indices)
        
        bool(pause_start_indices(j):pause_end_indices(j)) = false(1);
        
    end
    
    bool = bool(1:length(Dist_array));
    
    Dist_array = Dist_array(bool);
    
    velocities = Dist_array/(T*t0);
    
    velocities = velocities(velocities > Vmin_fit & velocities < Vmax_fit);
        
    dat_length(i) = length(velocities);

end

V_dat = zeros([L, max(dat_length)]);

for i=1:L
    
    filename = ['RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = selected_l0_array(i);
    
    factor = l0/l0_WLC;

    x_array = factor*trace_dat(:,8);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);
    
    len = length(xs_array);
    
    Dist_array = transpose(xs_array(T+1:len) - xs_array(1:len-T));
    
    [Tc_array] = find_crossing_time(t_array, xs_array,Nw);
    
    [DWT_array]=find_Dwell_time(Tc_array);
    
    DWT_indices = 1:length(DWT_array);
    
    pause_start_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V));
    
    pause_end_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V)+1);
    
    pause_start_indices = round(pause_start_times/t0 +1);
    
    pause_end_indices = round(pause_end_times/t0 +1);
    
    pause_start_indices = pause_start_indices-(T+1);
    
    pause_start_indices(pause_start_indices < 1) = 1;
    
    bool = true([1, len]);
    
    for j=1:length(pause_start_indices)
        
        bool(pause_start_indices(j):pause_end_indices(j)) = false(1);
        
    end
    
    bool = bool(1:length(Dist_array));
    
    Dist_array = Dist_array(bool);
    
    velocities = Dist_array/(T*t0);
    
    velocities = velocities(velocities > Vmin_fit & velocities < Vmax_fit);
        
    V_dat(i,1:dat_length(i)) = velocities;

end

save(str_V_dat, 'V_dat', '-mat')

save(str_V_length, 'dat_length', '-mat')
