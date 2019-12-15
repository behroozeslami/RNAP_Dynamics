W = 12;

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

nboot=100;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'Export_AF_7.5pN_1mM_NTP_LP3';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

DWT_folder = [output_foldername slash 'selected_DWT_pdf_Combined' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

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

str_DWT = [DWT_folder slash 'DWT_pdf_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

str_DWT_dat = [DWT_folder slash 'DWT_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

DWT_array = [];

L=length(selected_trace_number);

for j=1:L

    
    filename = ['RNAP_' num2str(selected_trace_number(j)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = selected_l0_array(j);
    
    factor = l0/l0_WLC;
    
    %factor = 1;

    x_array = factor*trace_dat(:,8);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);

    [Tc_array]=find_crossing_time(t_array, xs_array,Nw);

    [DWT_array_now] = find_Dwell_time(Tc_array);
        
    DWT_array = [DWT_array, DWT_array_now];

end

[bar_pos, bins]=make_bins(Ti,bins_per_decade,number_of_decades,t0,correction);

[DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);

[bins, bar_pos] = remove_empty_bins(bins, DWT_hist);

[DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);

bootstat = bootstrp(nboot,@ (DWT_array) DwellTimeHist_v3(DWT_array, t0, bins),DWT_array);

Mean_DWT_hist=zeros([1, length(DWT_hist)]);

Err=zeros([1, length(DWT_hist)]);

for i=1:length(DWT_hist)
    
    Mean_DWT_hist(:,i)=mean(bootstat(:,i));
    
    Err(:,i)=std(bootstat(:,i));
end

DWT_hist_dat = [bins', [bar_pos, NaN(1,1)]', [DWT_hist, NaN(1,1)]', [Err, NaN(1,1)]' ];

DWT_dat = DWT_array';

figure()

ax = axes();

PL = errorbar(bar_pos,DWT_hist,Err,'gO','MarkerSize',6, 'MarkerFaceColor', 'g');

set(ax,'XScale','log','YScale','log')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

xlabel('Dwell Time (s)')

ylabel('PDF')

saveas(PL, [str_DWT, '.jpg'], 'jpg')

saveas(PL, [str_DWT, '.fig'], 'fig')

save([str_DWT, '.txt'], 'DWT_hist_dat', '-ascii')

save([str_DWT_dat, '.txt'], 'DWT_dat', '-ascii')
