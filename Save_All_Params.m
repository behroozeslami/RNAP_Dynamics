W = 12;

k = 0;

Nw = 10;

Ts = 2*W+1;

T = Ts;

cam_freq = 25;

t0 = 1/cam_freq;

Pause_threshold = 100.64;

Pause_threshold_V = 100.64;

Dobcktrck = 0;

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

output_file_str = [output_foldername slash 'Final_Output' '_' output_foldername '.xlsx'];

folder_dat = transpose({'Folder name', output_foldername});

file_str_txt = [output_foldername slash 'Trace_data' '.txt'];

dat = load(file_str_txt, '-ascii');

trace_number = transpose(dat(:,1));

seidel_array = transpose(dat(:,2));

l0_array = transpose(dat(:,3));

l0_WLC = dat(1,4);

select_vec = (l0_array >= l0_min) & (l0_array <= l0_max) & (seidel_array >= Seidel_min) & (seidel_array <= Seidel_max);

selected_trace_number = floor(trace_number(select_vec));

dat = dat(select_vec,:);

Blacklist = load([output_foldername '/' 'Blacklist.txt'], '-ascii');

High_V_Blacklist = load([output_foldername '/' 'Blacklist_High_Velocity.txt'], '-ascii');

Low_V_Blacklist = load([output_foldername '/' 'Blacklist_Low_Velocity.txt'], '-ascii');

Blacklist = [Blacklist; High_V_Blacklist; Low_V_Blacklist];

index_2_keep = Exclude_Blacklist(selected_trace_number, Blacklist);

selected_trace_number = selected_trace_number(index_2_keep);

[sorted_selected_trace_number,Ind] = sort(selected_trace_number);

dat = dat(index_2_keep,:);

dat = dat(Ind,:);

DWT_pdf_selec_folder = [output_foldername slash 'DWTpdf_Selection' '_' subfolder '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

Max_Rel_dispersion = transpose(load([DWT_pdf_selec_folder '/' 'Max_Rel_Dispersion_in_kel.txt'], '-ascii'));

param_dat = [cam_freq, Ts*t0, T*t0, Nw, Pause_threshold, Pause_threshold_V, l0_min, l0_max, Max_Rel_dispersion];

Best_Fit_folder = [output_foldername slash 'Best_Fit' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_k_array = [Best_Fit_folder slash 'BestFit_k_array' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'T=' num2str(T) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_q_array = [Best_Fit_folder slash 'BestFit_q_array' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'T=' num2str(T) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

k_dat = load([str_k_array, '.txt'], '-ascii');

k_array = k_dat(1,:);

Low_Err_k = k_dat(2,:);

Up_Err_k = k_dat(3,:);

q_dat = load([str_q_array, '.txt'], '-ascii');

q_array = q_dat(1,:);

Low_Err_q = q_dat(2,:);

Up_Err_q = q_dat(3,:);

pause_num = length(q_array);

fit_dat = nan([pause_num+1, 8]);

fit_dat(:,1) = k_array;

fit_dat(:,2) = Low_Err_k;

fit_dat(:,3) = Up_Err_k;

fit_dat(1:pause_num,4) = q_array;

fit_dat(1:pause_num,5) = Low_Err_q;

fit_dat(1:pause_num,6) = Up_Err_q;

fit_dat(1,7) = round(pause_num);

fit_dat(1,8) = Dobcktrck;

header_1 = {'Trace number' 'Seidel (um)' 'EXP. conversion factor' 'WLC. conversion factor'}; 

header_2 = {'Camera frequency (Hz)' 'Filtering time (s)' 'Time window (s)' 'Dwell time window (bp)', 'Pause threshold for correction (s)', 'Pause threshold for exclusion (s)', 'Min. EXP. conversion factor (nm/bp)', 'Max. EXP conversion factor (nm/bp)', 'Max. dispersion in peak position (in unit sigma)'}; 

header_3 = {'Rates (Hz)', 'Low_Error in rates (Hz)', 'Up_Error in rates (Hz)', 'Pause probabilities', 'Low_Error in pause probabilities', 'Up_Error in pause probabilities', 'Number of Pauses', 'Backtrack included'};

xlswrite(output_file_str, folder_dat, 1)

xlswrite(output_file_str,dat,'Selected traces', 'A2')

xlswrite(output_file_str,header_1,'Selected traces', 'A1')

xlswrite(output_file_str,param_dat,'Parameters', 'A2')

xlswrite(output_file_str,header_2,'Parameters', 'A1');

xlswrite(output_file_str,fit_dat,'Best fit', 'A2')

xlswrite(output_file_str,header_3,'Best fit', 'A1');

table2word(header_2,param_dat, 'Table Grid 1', ['Parameters' '_' output_foldername '.docx'],'Parameters')

movefile(['Parameters' '_' output_foldername '.docx'], output_foldername)

fit_dat(isnan(fit_dat)) = 0;

table2word(header_3,fit_dat, 'Table Grid 1', ['BestFit' '_' output_foldername '.docx'],'Fit Results')

movefile(['BestFit' '_' output_foldername '.docx'], output_foldername)

