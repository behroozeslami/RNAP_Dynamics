W = 12;

k = 0;

Nw = 10;

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0 = 1/cam_freq;

Pause_threshold =  37.88;

t_strt = 0.1;

t_end = 10000;

nboot = 100;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'Select_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

%DWT_folder = ['Modelling' slash output_foldername];
DWT_folder = output_foldername;

if ~(exist(DWT_folder,'dir')==7)
    
    mkdir(DWT_folder)
    
end

DWT_folder_out = [output_foldername slash 'selected_DWT_pdf_Corrected_Single_Traces' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

if ~(exist(DWT_folder_out,'dir')==7)
    
    mkdir(DWT_folder_out)
    
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

str_DWT_length = [DWT_folder slash 'DWT_data_length_Corrected' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's' '.mat'];

str_DWT_dat = [DWT_folder slash 'DWT_data_Corrected' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's' '.mat'];

str_DWT_bins = [DWT_folder slash 'DWT_pdf_bincenters_Corrected' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's' '.mat'];

str_DWT_hist = [DWT_folder slash 'DWT_pdf_Corrected' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's' '.mat'];

str_DWT_error = [DWT_folder slash 'DWT_pdf_error_Corrected' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's' '.mat'];

str_trace_index = [DWT_folder slash 'Trace_indices' '.mat'];

DWT_folder_0 = [output_foldername slash 'selected_DWT_pdf_Corrected' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_DWT_0 = [DWT_folder_0 slash 'DWT_pdf_Corrected' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

DWT_hist_dat_0 = load([str_DWT_0, '.txt'], '-ascii');

bins = transpose(DWT_hist_dat_0(:,1));

centers = transpose(DWT_hist_dat_0(1:length(bins)-1,2));

DWT_hist0 = transpose(DWT_hist_dat_0(1:length(bins)-1,3));

Err0 = transpose(DWT_hist_dat_0(1:length(bins)-1,4));

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
    
    Tc_array = correct_crossing_time(Tc_array, Pause_threshold, t0);

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
    
    Tc_array = correct_crossing_time(Tc_array, Pause_threshold, t0);

    [DWT_array] = find_Dwell_time(Tc_array);
    
    bootstat = bootstrp(nboot,@ (DWT_array) DwellTimeHist_v3(DWT_array, t0, bins),DWT_array);
    
    Mean_DWT_hist = mean(bootstat);
    
    Err = std(bootstat);
    
    DWT_hist_dat(i,:) = Mean_DWT_hist;
    
    DWT_hist_err(i,:) = Err;
    
    DWT_array = DWT_array(DWT_array > t_strt & DWT_array < t_end);
        
    DWT_dat(i,1:dat_length(i)) = DWT_array;
    
    str_DWT_hist_out = [DWT_folder_out slash 'Trace' '_' num2str(selected_trace_number(i)) '_' 'DWT_pdf_Corrected' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];
    
    figure()

    ax = axes();

    PL = errorbar(centers,DWT_hist0,Err0,'kO','MarkerSize',6, 'MarkerFaceColor', 'k');
    
    hold on
    
    errorbar(centers,Mean_DWT_hist,Err,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');
    
    hold off

    set(ax,'XScale','log','YScale','log')

    str_title = ['Trace' num2str(selected_trace_number(i)) ' ' 'Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw)];

    title(str_title)

    xlabel('Dwell Time (s)')

    ylabel('PDF')

    saveas(PL, [str_DWT_hist_out, '.jpg'], 'jpg')

    saveas(PL, [str_DWT_hist_out, '.fig'], 'fig')
    
    close()

end

save(str_DWT_dat, 'DWT_dat', '-mat')

save(str_DWT_length, 'dat_length', '-mat')

save(str_DWT_hist, 'DWT_hist_dat', '-mat')

save(str_DWT_error, 'DWT_hist_err', '-mat')

save(str_trace_index, 'selected_trace_number', '-mat');

save(str_DWT_bins, 'centers', '-mat')
