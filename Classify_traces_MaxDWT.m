W = 24;

k = 0;

Nw = 4;

cam_freq = 25;
%%%%%% IMPORTANT
cam_freq = 50;

t0 = 1/cam_freq;

t_strt = 0.05;

t_end = 10000;

mainfoldername = 'OF4_9pN_1mM_NTP_GreA';

data_name = 'GreA';

compare_with3 = logical(1);

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

directory = ['..' slash 'Richard_Data' slash mainfoldername slash subfolder];

if strcmp('',subfolder)
    
    directory = ['..' slash 'Richard_Data' slash mainfoldername];
    
end

str_trace_index = [DWT_folder slash 'Trace_indices' '.mat'];

trace_indices = load(str_trace_index);

trace_indices = trace_indices.selected_trace_number;

str_DWT_length = [DWT_folder slash 'DWT_data_length' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

dat_length = load(str_DWT_length);

DWT_length = dat_length.dat_length;

str_DWT_dat = [DWT_folder slash 'DWT_data' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

DWTs = load(str_DWT_dat);

DWTs = DWTs.DWT_dat;

Max_DWT = arrayfun(@(row) max(DWTs(row, :)), 1:length(trace_indices));

opt = statset('MaxIter',1000);

GMM_2 = fitgmdist(Max_DWT',2,'Options', opt);

GMM_3 = fitgmdist(Max_DWT',3,'Options', opt);

GMM = GMM_2;

G_num = 2;

names ={'S_DWT', 'L_DWT'};

if GMM_3.BIC < GMM_2.BIC && compare_with3
    
    GMM = GMM_3;
    
    G_num = 3;
    
    names ={'S_DWT', 'M_DWT', 'L_DWT'};
    
end

mu = GMM.mu;

[sorted_mu, I] = sort(mu);

ids = 1:G_num;

ids = ids(I);

idx = cluster(GMM, Max_DWT'); 

select_folder = [DWT_folder slash 'Selection_maxDWT_GMM' num2str(G_num)];

if ~(exist(select_folder,'dir')==7)
    
    mkdir(select_folder)
    
end

for n=1:G_num

    selected_trace_indices = trace_indices(idx==ids(n));

    name = names{n};

    class_name = [select_folder slash data_name, '_', name '.txt'];

    save(class_name, 'selected_trace_indices', '-ascii')

end

if G_num == 3

    selected_trace_indices = trace_indices(~(idx==ids(G_num)));

    name = 'SM_DWT';

    class_name = [select_folder slash data_name, '_', name '.txt'];

    save(class_name, 'selected_trace_indices', '-ascii')
    
end


Max_DWT_bins = 0:10:ceil(10*max(Max_DWT))/10;

colors = {'g','r','b'};

if G_num ==2

    colors = {'r','b'};

end

for n=1:G_num

PL = histogram(Max_DWT(idx==ids(n)), Max_DWT_bins, 'EdgeColor','k', 'FaceColor', colors{n});

hold on

end

hold off

xlabel('Max. Dwell Time (s)')

ylabel('count')

title(data_name)

str_fig = [select_folder slash data_name, '_' 'MaxDWT_selection'];

saveas(PL, [str_fig, '.bmp'], 'bmp')

saveas(PL, [str_fig, '.fig'], 'fig')

str_model = [select_folder slash 'GMM' num2str(G_num) '_' data_name];

save([str_model, '.mat'], 'GMM', '-mat')

disp([num2str(G_num) ' Gaussians were fitted.'])