W = 12;

k = 0;

Nw = 10;

cam_freq = 25;

t0 = 1/cam_freq;

%Pause_threshold =  37.88;

Pause_threshold =  10;

T = round(Pause_threshold/t0);

bins_per_decade=7;

Ti=0.01;

number_of_decades=5;

correction=1;

nboot=100;

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

trace_folder = [output_foldername slash 'velocity_traces' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

if ~(exist(trace_folder,'dir')==7)
    
    mkdir(trace_folder)
    
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

figure()

L = length(selected_trace_number);

%L = 1;

for i=1:L
    
    filename = ['RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    t_array = t_array - t_array(1);
    
    l0 = selected_l0_array(i);
    
    factor = l0/l0_WLC;

    x_array = factor*trace_dat(:,8);
    
    [xs_array] = Smooth_Trace_SG(x_array,W, k);
    
    len = length(xs_array);
    
    Dists = transpose(xs_array(T+1:len) - xs_array(1:len-T));
    
    velocities = Dists/(T*t0);
    
    new_t_array = (0:length(velocities)-1)*t0 + 0.5*T*t0;
    
    PL = plot(new_t_array, velocities, 'g', 'LineWidth', 2);
    
    hold on
    
    [Tc_array] = find_crossing_time(t_array, xs_array,Nw);
    
    [DWT_array]=find_Dwell_time(Tc_array);
    
    DWT_indices = 1:length(DWT_array);
    
    pause_start_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold));
    
    pause_end_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold)+1);
    
    pause_start_indices = round(pause_start_times/t0 +1);
    
    pause_end_indices = round(pause_end_times/t0 +1);
    
    pause_bool_vec = ones([1, length(t_array)]);
    
    for n=1:length(pause_start_indices)
        
        pause_start = pause_start_indices(n);
        
        pause_end = pause_end_indices(n);
        
        pause_bool_vec(pause_start:pause_end) = 0;
    end
    
    delta_pause_bool_vec = pause_bool_vec(2:length(pause_bool_vec)) - pause_bool_vec(1:length(pause_bool_vec)-1);
    
    indices = 1:length(t_array);
    
    chunk_start = indices(delta_pause_bool_vec == 1);
    
    chunk_end = indices(delta_pause_bool_vec == -1);
    
    if (length(chunk_start) > 0) && (length(chunk_end) > 0)  
        
        if chunk_start(1) > chunk_end(1)
        
            chunk_start =  [1, chunk_start];
        end
        
        if chunk_start(length(chunk_start)) > chunk_end(length(chunk_end))
        
            chunk_end =  [chunk_end, indices(length(indices))];
        end
        
    else
        
        chunk_start = 1;
        
        chunk_end = length(t_array);
        
    end
    
    remove_end_chunks = chunk_start<=length(velocities);
    
    chunk_end = chunk_end(remove_end_chunks);
    
    chunk_start = chunk_start(remove_end_chunks);
    
    if chunk_end(length(chunk_end)) > length(velocities)
        
        chunk_end(length(chunk_end)) = length(velocities);
        
    end
    
    for j=1:length(chunk_start)
        
        v_chunk = velocities(chunk_start(j):chunk_end(j));
        
        t_chunk = t_array(chunk_start(j):chunk_end(j))+0.5*T*t0;
       
        
        plot(t_chunk, v_chunk, 'r', 'LineWidth', 2)
    end
    
    for j=1:length(pause_start_times)
        
        line([pause_start_times(j), pause_start_times(j)], ylim)
    end
    
    for j=1:length(pause_end_times)
        
        line([pause_end_times(j), pause_end_times(j)], ylim)
    end
    
    hold off
    
    xlabel('Time (s)')

    ylabel('Position (bp)')
    
    str_trace =[trace_folder slash 'RNAP_' num2str(selected_trace_number(i))];
    
    saveas(PL, [str_trace, '.jpg'], 'jpg')

    %saveas(PL, [str_trace, '.fig'], 'fig')

end

