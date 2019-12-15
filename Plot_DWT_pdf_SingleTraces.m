W = 12;

k = 0;

Nw = 4;

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0 = 1/cam_freq;

bins_per_decade=7;

Ti=0.01;

number_of_decades=6;

correction=1;

nboot=100;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

%mainfoldername = 'Export_AF_7.5pN_1mM_NTP_I';
mainfoldername = 'EXPORT_OF_12.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

DWT_folder = [output_foldername slash 'selected_DWT_pdf_Single_Traces' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

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

str_DWT = [DWT_folder slash 'DWT_pdf_Single_Traces' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

str_DWTmax = [DWT_folder slash 'MaxDWT_hist_Single_Traces' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

L = length(selected_trace_number);

C_L = 3;

cmap = [0,0,1;0,1,0;1,0,0];

DWT_max_array = zeros([1,L]);

Group_1 = [];

Group_2 = [];

Group_3 = [];

PL = figure();

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
    
    [bar_pos, bins]=make_bins(Ti,bins_per_decade,number_of_decades,t0,correction);
    
    [DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);

    %[bins, bar_pos] = remove_empty_bins(bins, DWT_hist);

    %[DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);
    
    %DWT_max = max(bar_pos);
    
    DWT_max = max(DWT_array);
    
    DWT_max_array(i) = DWT_max;
    
    C_index = ceil(log10(DWT_max));
    
    C_index = max([C_index, 1]);
    
    C_index = min([C_index, C_L]);
    
    if C_index == 1
        
        Group_1 = [Group_1, selected_trace_number(i)];
        
    end
    
    if C_index == 2
        
        Group_2 = [Group_2, selected_trace_number(i)];
        
    end
    
    if C_index == 3
        
        Group_3 = [Group_3, selected_trace_number(i)];
        
    end
    
    loglog(bar_pos,DWT_hist,'O-','MarkerSize',6, 'Color', cmap(C_index,:), 'MarkerFaceColor', cmap(C_index,:));
    
    hold on

end

hold off

str_title = ['Ts = ' ' ' num2str((2*W+1)*t0) 's' ' ' 'Nw = ' num2str(Nw)];

title(str_title)

xlabel('Dwell Time (s)')

ylabel('PDF')

saveas(PL, [str_DWT, '.jpg'], 'jpg')

saveas(PL, [str_DWT, '.fig'], 'fig')

PL2 = figure();

centers = -2:4;

h = hist(log10(DWT_max_array), centers);

bar(centers, h, 'Facecolor', 'r');

xlabel('Log10(Max Dwell Time)')

ylabel('number of traces')

saveas(PL2, [str_DWTmax, '.jpg'], 'jpg')

saveas(PL2, [str_DWTmax, '.fig'], 'fig')

save([output_foldername '/' 'TraceID_ShortPause.txt'], 'Group_1', '-ascii')

save([output_foldername '/' 'TraceID_IntermediatePause.txt'], 'Group_2', '-ascii')

save([output_foldername '/' 'TraceID_LongPause.txt'], 'Group_3', '-ascii')

