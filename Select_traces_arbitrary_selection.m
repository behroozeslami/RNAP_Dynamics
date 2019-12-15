W = 25;

Ts = 2*W + 1;

k = 0;

Nw = 4;

%cam_freq = 25;
%%%%%% IMPORTANT
cam_freq = 50;

t0 = 1/cam_freq;

bins_per_decade=7;

Ti=0.01;

number_of_decades=6;

correction=1;

Ti_peak = 0.1;

Tf_peak = Ts*t0;

nboot = 1000;

min_App_kel = 8;

max_App_kel = 15;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'SELECTED_OF4_9pN_1mM_NTP_GreA_2';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

DWT_pdf_selec_folder = [output_foldername slash 'DWTpdf_Selection' '_' subfolder '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

if ~(exist(DWT_pdf_selec_folder,'dir')==7)
    mkdir(DWT_pdf_selec_folder)
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

[bar_pos, bins]=make_bins(Ti,bins_per_decade,number_of_decades,t0,correction);

DWT_array = [];

App_kel = zeros([length(selected_trace_number), 1]);

Err_App_kel = zeros([length(selected_trace_number), 1]);

min_DWT = Inf;

max_DWT = 0;

for i=1:length(selected_trace_number)
    
    filename = ['RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = selected_l0_array(i);    
    
    factor = l0/l0_WLC;

    x_array = factor*trace_dat(:,8);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);

    [Tc_array]=find_crossing_time(t_array, xs_array,Nw);

    [DWT_array]=find_Dwell_time(Tc_array);
    
    out = get_DWTpdf_params(DWT_array, Ti_peak, Tf_peak);
    
    App_kel(i,1) = Nw/out(1);
    
    bootstat_param = bootstrp(nboot,@ (DWT_array) get_DWTpdf_params(DWT_array, Ti_peak, Tf_peak), DWT_array);

    App_k_el_array = Nw./bootstat_param(:,1);
    
    Err_App_kel(i,1) = sqrt(mean((App_k_el_array - App_kel(i,1)).^2));
    
    min_DWT_now = min(DWT_array);
    
    min_DWT = min([min_DWT, min_DWT_now]);
    
    max_DWT_now = max(DWT_array);
    
    max_DWT = max([max_DWT, max_DWT_now]);

end


Low_V_Blacklist = transpose(selected_trace_number((App_kel < min_App_kel) | isnan(App_kel)));

High_V_Blacklist = transpose(selected_trace_number((App_kel > max_App_kel)));

save([output_foldername '/' 'Blacklist_High_Velocity.txt'], 'High_V_Blacklist', '-ascii')

save([output_foldername '/' 'Blacklist_Low_Velocity.txt'], 'Low_V_Blacklist', '-ascii')


