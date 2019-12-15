cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0 = 1/cam_freq;

bins_per_decade=7;

Ti=0.01;

number_of_decades=6;

correction=1;

nboot=100;

mainfoldername = 'NEW_AF_7.5pN_1mM_ITP_100uM_NTP';

file_name = 'Control_dissociation';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

DWT_folder = [output_foldername slash 'EndPause_pdf'];

if ~(exist(DWT_folder,'dir')==7)
    
    mkdir(DWT_folder)
    
end

DWT_array = xlsread([output_foldername slash file_name '.xlsx']);

[bar_pos, bins]=make_bins(Ti,bins_per_decade,number_of_decades,t0,correction);

[DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);

%[bins, bar_pos] = remove_empty_bins(bins, DWT_hist);

[DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);

bootstat = bootstrp(nboot,@ (DWT_array) DwellTimeHist_v3(DWT_array, t0, bins),DWT_array);

Mean_DWT_hist=zeros([1, length(DWT_hist)]);

Err=zeros([1, length(DWT_hist)]);

for i=1:length(DWT_hist)
    
    Mean_DWT_hist(:,i)=mean(bootstat(:,i));
    
    Err(:,i)=std(bootstat(:,i));
end

figure()

ax = axes();

PL = errorbar(bar_pos,DWT_hist,Err,'kO','MarkerSize',6, 'MarkerFaceColor', 'k');

set(ax,'XScale','log','YScale','log');

title(strrep(file_name,'_', '  '))

xlabel('Dwell Time (s)')

ylabel('PDF')

str_DWT = [DWT_folder slash file_name];

saveas(PL, [str_DWT, '.jpg'], 'jpg')

saveas(PL, [str_DWT, '.fig'], 'fig')

