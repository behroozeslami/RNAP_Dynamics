W = 12;

k = 0;

Nw = 10;

cam_freq = 25;

t0 = 1/cam_freq;

t_strt = 0.05;

t_end = 10000;

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

DWT_folder = ['Modelling' slash output_foldername];

if ~(exist(DWT_folder,'dir')==7)
    
    mkdir(DWT_folder)
    
end

str_DWT_dat = [DWT_folder slash 'DWT_data' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

str_DWT_bincenters = [DWT_folder slash 'DWT_pdf_bincenters' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];

str_DWT_bins = [DWT_folder slash 'DWT_pdf_bins' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];

str_DWT_hist = [DWT_folder slash 'DWT_pdf' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];


DWT_dat = load(str_DWT_dat, '-mat');

DWT_dat = DWT_dat.DWT_dat';

DWT_hist_dat = load(str_DWT_hist, '-mat');

DWT_hist_dat = DWT_hist_dat.DWT_hist_dat';

centers = load(str_DWT_bincenters, '-mat');

centers = centers.centers;

bins = load(str_DWT_bins, '', '-mat');

bins = bins.bins;

max_DWT = max(DWT_dat);

H_folder = [DWT_folder slash 'Heterogeneity'];

if ~(exist(H_folder,'dir')==7)
    
    mkdir(H_folder)
    
end

figure()

ax = axes();

% long DWT_pdf

data = [centers', DWT_hist_dat(:, max_DWT>100)];

colnames = {'bin centers'};

for j=2:length(data)
    
    colnames{j} = ['hist_' num2str(j-1)];
    
    loglog(centers, data(:,j), 'r.', 'MarkerSize', 8)
    
    hold on
    
end

str_DWT_out = [H_folder slash 'Long_DWT_pdf' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.xlsx'];

xlswrite(str_DWT_out,data,'Sheet1','A2');

xlswrite(str_DWT_out,colnames,'Sheet1','A1');

% Medium DWT_pdf

data = [centers', DWT_hist_dat(:, max_DWT>10 &max_DWT<100)];

colnames = {'bin centers'};

for j=2:length(data)
    
    colnames{j} = ['hist_' num2str(j-1)];
    
    loglog(centers, data(:,j), 'g.', 'MarkerSize', 8)
    
    hold on
    
end

str_DWT_out = [H_folder slash 'Medium_DWT_pdf' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.xlsx'];

xlswrite(str_DWT_out,data,'Sheet1','A2');

xlswrite(str_DWT_out,colnames,'Sheet1','A1');

% Short DWT_pdf

data = [centers', DWT_hist_dat(:, max_DWT<10)];

colnames = {'bin centers'};

for j=2:length(data)
    
    colnames{j} = ['hist_' num2str(j-1)];
    
    loglog(centers, data(:,j), 'b.', 'MarkerSize', 8)
    
    hold on
    
end

str_DWT_out = [H_folder slash 'Short_DWT_pdf' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.xlsx'];

xlswrite(str_DWT_out,data,'Sheet1','A2');

xlswrite(str_DWT_out,colnames,'Sheet1','A1');

DWT_folder_0 = [output_foldername slash 'selected_DWT_pdf_Combined' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

str_DWT_0 = [DWT_folder_0 slash 'DWT_pdf_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

DWT_hist_dat_0 = load([str_DWT_0, '.txt'], '-ascii');

DWT_hist = transpose(DWT_hist_dat_0(1:length(bins)-1,3));

Err = transpose(DWT_hist_dat_0(1:length(bins)-1,4));

PL = errorbar(centers,DWT_hist,Err,'kO','MarkerSize',6, 'MarkerFaceColor', 'k');

set(ax,'XScale','log','YScale','log')

str_title = ['Ts = ' ' ' num2str((2*W+1)*t0) 's' ' ' 'Nw = ' num2str(Nw)];

title(str_title)

xlabel('Dwell Time (s)')

ylabel('PDF')

hold off

str_DWT_fig = [H_folder slash 'DWT_pdf' '_' '1mM7AF' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

saveas(PL, [str_DWT_fig, '.jpg'], 'jpg')

saveas(PL, [str_DWT_fig, '.fig'], 'fig')
