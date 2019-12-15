BIC_shift = 1;

f0 = 12;

f = 7.5;

B = exp(-f/f0);

W = 12;

k = 0;

Nw = 10;

cam_freq = 25;

t0 = 1/cam_freq;

Pause_threshold =  100.64;

bins_per_decade = 7;

Ti = 0.01;

number_of_decades = 6;

correction = 0;

t_strt = 0.1;

t_end = 10000;

pause_num_array = [2, 3, 4, 5];

fitCount_tot = 5;

BIC_dat = zeros([length(pause_num_array), 4]);

Dobcktrck = 0;

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

Fit_folder = [output_foldername slash 'BIC_Results' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

if ~(exist(Fit_folder,'dir')==7)
    
    mkdir(Fit_folder)
    
end

best_fit_num = zeros([1, length(pause_num_array)]);

for n=1:length(pause_num_array)
    
    pause_num = pause_num_array(n);
    
    Min_ML = Inf;
    
    Blacklist = load([Fit_folder slash 'Blacklist_Pause_num=' num2str(pause_num) '.txt'], '-ascii');
    
    v = true([1, fitCount_tot]);
    
    v(Blacklist) = false;
    
    fit_nums = 1:fitCount_tot;
    
    fit_nums = fit_nums(v);
    
    for j=1:length(fit_nums)
        
        fit_num = fit_nums(j);
        
        str_fit = [Fit_folder slash 'Fit_Result' '_' 'Pause_num=' num2str(pause_num) '_' 'fit_num=' num2str(fit_num) '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];
        
         fit_dat = load([str_fit, '.txt'], '-ascii');
         
         fval = fit_dat(1,3);
         
         if Min_ML > fval
            
            Min_ML = fval;
            
            BIC = fit_dat(2,3);
            
            best_fit_num(n) = fit_num;
            
        end
    end
    
    BIC_dat(n, :) = [pause_num, best_fit_num(n), Min_ML, BIC];
end

str_BIC = [Fit_folder slash 'BIC_Result' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

save([str_BIC, '.txt'], 'BIC_dat', '-ascii')

BIC_array = transpose(BIC_dat(:,4));

delta_BIC_array = BIC_array - min(BIC_array) + BIC_shift;

figure()

PL = semilogy(pause_num_array, delta_BIC_array, 'rO-','MarkerSize',6, 'MarkerFaceColor', 'r')

xlabel('Number of Pauses')

ylabel('\Delta BIC')

str_title = ['BT = ' num2str(Dobcktrck) ' ' 'Ts = ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];
    
title(str_title)

set(gca,'XTick',pause_num_array, 'YMinorTick','on')

grid on

grid minor

saveas(PL, [str_BIC, '.jpg'], 'jpg')

saveas(PL, [str_BIC, '.fig'], 'fig')
