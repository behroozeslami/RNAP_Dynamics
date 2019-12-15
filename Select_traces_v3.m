W = 12;

Ts = 2*W + 1;

k = 0;

Nw = 10;

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0 = 1/cam_freq;

bins_per_decade=7;

Ti=0.01;

number_of_decades=6;

correction=1;

Ti_peak = 0.1;

Tf_peak = Ts*t0;

Max_Rel_dispersion = 2;

nboot = 1000;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'Export_AF_7.5pN_1mM_NTP_Long';

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

max_DWT = min([max_DWT, max(bins)]);

mean_Err = nanmean(Err_App_kel);

binwidth = mean_Err;

centers = min(App_kel)-binwidth:binwidth:max(App_kel)+binwidth;

App_kel_hist = hist(App_kel, centers);

[Low_Bond, High_Bond] = Tukey(App_kel, 1.5);

best_App_kel = nanmean(App_kel(App_kel > Low_Bond & App_kel < High_Bond));

[min_dist, peak_index] = min(abs(centers- best_App_kel));

std_App_kel = nanstd(App_kel(App_kel > Low_Bond & App_kel < High_Bond));

selection_param =  abs(App_kel - best_App_kel)./(Max_Rel_dispersion*std_App_kel);

selected_App_kel = App_kel(selection_param <= 1);

select_App_kel_hist = hist(selected_App_kel, centers);

Low_V_Blacklist = transpose(selected_trace_number((selection_param > 1 & App_kel < best_App_kel) | isnan(App_kel)));

High_V_Blacklist = transpose(selected_trace_number((selection_param > 1 & App_kel > best_App_kel)));

figure()

PL3 = bar(centers, App_kel_hist, 'Facecolor', 'b', 'Edgecolor', 'k');

hold on

bar(centers, select_App_kel_hist, 'Facecolor', 'r', 'Edgecolor', 'k');

line([Low_Bond, Low_Bond], ylim, 'Color','g','LineWidth',2)

line([High_Bond, High_Bond], ylim, 'Color','g','LineWidth',2)

hold off

xlabel('Apparent elongation rate (bp/s)')

ylabel('Counts (number of traces)')

legend('All', 'Selected')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

str_kel = [DWT_pdf_selec_folder slash 'Apparent_kel_hist' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

saveas(PL3, [str_kel, '.jpg'], 'jpg')

saveas(PL3, [str_kel, '.fig'], 'fig')

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
    
    markersize = 3;
    
    linewidth = 1;
    
    if (selection_param(i,1) <= 1)
        
        markersize = 6;
    
        linewidth = 3;
    
    end
    
    index = indices(centers <= (App_kel_now + 0.5*binwidth) & centers > (App_kel_now - 0.5*binwidth));
    
    c_index = I-(abs(index - peak_index)); 
    
    if (isnan(App_kel_now))
        
        c_index = 1;
        
        markersize = 3;
    
        linewidth = 1;
        
    end
    
    [DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);
    
    PL4 = loglog(bar_pos,DWT_hist,'O-','LineWidth',linewidth, 'Color', cmap(c_index,:), 'MarkerSize',markersize, 'MarkerFaceColor', cmap(c_index,:));
    
    hold on

end

hold off

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

xlim([min_DWT, max_DWT])

xlabel('Dwell Time (s)')

ylabel('PDF')

str_DWT_pdf = [DWT_pdf_selec_folder slash 'DWT_pdf_Selected' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

saveas(PL4, [str_DWT_pdf, '.jpg'], 'jpg')

saveas(PL4, [str_DWT_pdf, '.fig'], 'fig')

save([output_foldername '/' 'Blacklist_High_Velocity.txt'], 'High_V_Blacklist', '-ascii')

save([output_foldername '/' 'Blacklist_Low_Velocity.txt'], 'Low_V_Blacklist', '-ascii')

save([DWT_pdf_selec_folder '/' 'Max_Rel_Dispersion_in_kel.txt'], 'Max_Rel_dispersion', '-ascii')

