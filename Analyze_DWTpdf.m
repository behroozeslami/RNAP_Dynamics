W = 12;

Ts = 2*W + 1;

k = 0;

Nw = 10;

cam_freq = 25;

t0 = 1/cam_freq;

bins_per_decade=7;

Ti=0.01;

number_of_decades=5;

correction=1;

Ti_peak = 0.1;

Tf_peak = Ts*t0;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'EXPORT_AF_7.5pN_0.5mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

DWT_pdf_param_folder = [output_foldername slash 'DWTpdf_parameters' '_' subfolder '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

if ~(exist(DWT_pdf_param_folder,'dir')==7)
    mkdir(DWT_pdf_param_folder)
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

selected_seidel_array = seidel_array(select_vec);

[bar_pos, bins]=make_bins(Ti,bins_per_decade,number_of_decades,t0,correction);

DWT_array = [];

DWTpdf_params = zeros([length(selected_trace_number), 3]);

sigma_array = zeros([1, length(selected_trace_number)]);

min_DWT = Inf;

max_DWT = 0;

for i=1:length(selected_trace_number)
    
    filename = ['RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = selected_l0_array(i);    
    
    factor = l0/l0_WLC;
    
    %factor = 1;

    x_array = factor*trace_dat(:,8);
    
    sigma_array(i) = Estimate_fast_noise(x_array);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);

    [Tc_array]=find_crossing_time(t_array, xs_array,Nw);

    [DWT_array]=find_Dwell_time(Tc_array);
    
    out = get_DWTpdf_params(DWT_array, Ti_peak, Tf_peak);
    
    DWTpdf_params(i,:) = out;
    
    min_DWT_now = min(DWT_array);
    
    min_DWT = min([min_DWT, min_DWT_now]);
    
    max_DWT_now = max(DWT_array);
    
    max_DWT = max([max_DWT, max_DWT_now]);

end

max_DWT = min([max_DWT, max(bins)]);

peak_pos = DWTpdf_params(:,1);

sigma_LG = DWTpdf_params(:,2);

P = DWTpdf_params(:,3);

App_kel = Nw ./peak_pos;

binwidth = 1;

centers = min(App_kel)-binwidth:binwidth:max(App_kel)+binwidth;

App_kel_hist = hist(App_kel, centers);

figure()

PL1 = plot(peak_pos,P, 'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

xlabel('Peak position (s)')

ylabel('Cumulative probability')

str_peak1 = [DWT_pdf_param_folder slash 'peak_vs_pauseprob' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

saveas(PL1, [str_peak1, '.jpg'], 'jpg')

saveas(PL1, [str_peak1, '.fig'], 'fig')

figure()

PL2 = plot(peak_pos,sigma_LG, 'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

xlabel('Peak position (s)')

ylabel('Lognormal Std (s)')

str_peak2 = [DWT_pdf_param_folder slash 'peak_vs_sigma' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

saveas(PL2, [str_peak2, '.jpg'], 'jpg')

saveas(PL2, [str_peak2, '.fig'], 'fig')

figure()

PL3 = bar(centers, App_kel_hist, 'Facecolor', 'r', 'Edgecolor', 'k');

xlabel('Apparent elongation rate (bp/s)')

ylabel('Counts (number of traces)')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

str_kel = [DWT_pdf_param_folder slash 'Apparent_kel_hist' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

saveas(PL3, [str_kel, '.jpg'], 'jpg')

saveas(PL3, [str_kel, '.fig'], 'fig')

figure()

PL5 = plot(App_kel,sigma_array, 'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

xlabel('Apparent elongation rate (bp/s)')

ylabel('Noise Amplitude (bp)')

str_sigma = [DWT_pdf_param_folder slash 'Apparent_kel_vs_sigma' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

saveas(PL5, [str_sigma, '.jpg'], 'jpg')

saveas(PL5, [str_sigma, '.fig'], 'fig')

figure()

PL6 = plot(App_kel,selected_l0_array, 'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

xlabel('Apparent elongation rate (bp/s)')

ylabel('Exp. conv. factor (nm/bp)')

str_l0 = [DWT_pdf_param_folder slash 'Apparent_kel_vs_Exp. conv. factor' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

saveas(PL6, [str_l0, '.jpg'], 'jpg')

saveas(PL6, [str_l0, '.fig'], 'fig')

figure()

PL7 = plot(App_kel,selected_seidel_array, 'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

xlabel('Apparent elongation rate (bp/s)')

ylabel('seidel (\mum)')

str_seidel = [DWT_pdf_param_folder slash 'Apparent_kel_vs_seidel' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

saveas(PL7, [str_seidel, '.jpg'], 'jpg')

saveas(PL7, [str_seidel, '.fig'], 'fig')

[max_hist, peak_index] = max(App_kel_hist);

best_App_kel = centers(peak_index);

I = max([peak_index - 1, length(centers) - peak_index]) + 1;

cmap = jet(I);

indices = 1:length(centers);

figure()

for i=1:length(selected_trace_number)
    
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

    [DWT_array]=find_Dwell_time(Tc_array);
    
    out = get_DWTpdf_params(DWT_array, Ti_peak, Tf_peak);
    
    App_kel_now = Nw/out(1);
    
    index = indices(centers <= (App_kel_now + 0.5*binwidth) & centers > (App_kel_now - 0.5*binwidth));
    
    c_index = I-(abs(index - peak_index)); 
    
    if (isnan(App_kel_now))
        
        c_index = 1;
        
    end
        
    [DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);
    
    PL4 = loglog(bar_pos,DWT_hist,'O-','LineWidth',3, 'Color', cmap(c_index,:), 'MarkerSize',6, 'MarkerFaceColor', cmap(c_index,:));
    
    hold on

end

hold off

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

xlim([min_DWT, max_DWT])

xlabel('Dwell Time (s)')

ylabel('PDF')

str_DWT_pdf = [DWT_pdf_param_folder slash 'DWT_pdf' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

saveas(PL4, [str_DWT_pdf, '.jpg'], 'jpg')

saveas(PL4, [str_DWT_pdf, '.fig'], 'fig')



