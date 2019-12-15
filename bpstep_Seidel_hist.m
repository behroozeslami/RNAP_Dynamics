
mainfoldername = 'Export_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end


if ~(exist(output_foldername,'dir')==7)
    
    mkdir(output_foldername)
    
end

hist_folder = [output_foldername slash 'Histograms_and_Correlations'];

if ~(exist(hist_folder,'dir')==7)
    
    mkdir(hist_folder)
    
end

directory = ['..' slash 'Richard_Data' slash mainfoldername slash subfolder];

if strcmp('',subfolder)
    directory = ['..' slash 'Richard_Data' slash mainfoldername];
end

file_str_txt = [output_foldername slash 'Trace_data' '.txt'];

dat = load(file_str_txt, '-ascii');

l0_array = transpose(dat(:,3));

l0_WLC = dat(1,4);

seidel_array = transpose(dat(:,2)); 

% bp_step_hist

figure()

ave_l0 = mean(l0_array);

std_l0 = std(l0_array);

binwidth = 0.01;

centers = min(l0_array):binwidth:max(l0_array);

h = hist(l0_array, centers);

h = h/(binwidth*sum(h));

PL1 = bar(centers, h, 'Facecolor', 'r', 'Edgecolor', 'k');

hold on

line([l0_WLC, l0_WLC], ylim, 'LineStyle', '--', 'LineWidth',2, 'Color', 'b')

hold off

xlabel('base pair step (nm)')

ylabel('PDF')

str_bps = [hist_folder slash 'bps_PDF'];

saveas(PL1, str_bps, 'jpg')
saveas(PL1, str_bps, 'fig')

% Seidel hist

figure()

binwidth = 0.01;

centers = min(seidel_array):binwidth:max(seidel_array);

h = hist(seidel_array, centers);

h = h/(binwidth*sum(h));

PL3 = bar(centers, h, 'Facecolor', 'r', 'Edgecolor', 'k');

xlabel('Seidel (\mu m)')

ylabel('PDF')

str_seidel = [hist_folder slash 'seidel_PDF'];

saveas(PL3, str_seidel, 'jpg')

saveas(PL3, str_seidel, 'fig')

% Correlations

figure()

PL4 = plot(seidel_array, l0_array, 'r.', 'MarkerSize', 12);

xlabel('Seidel (\mu m)')

ylabel('bp step (bp)')

str_seidel_bpstep = [hist_folder slash 'seidel_bpstep'];

saveas(PL4, str_seidel_bpstep, 'jpg')

saveas(PL4, str_seidel_bpstep, 'fig')





