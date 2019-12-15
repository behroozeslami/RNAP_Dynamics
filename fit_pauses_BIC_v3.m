f0 = 12;

f = 7.5;

B = exp(-f/f0);

W = 12;

k = 0;

Nw = 10;

cam_freq = 25;

t0 = 1/cam_freq;

%Pause_threshold =  37.88;
Pause_threshold =  100.64;

bins_per_decade = 7;

Ti = 0.01;

number_of_decades = 6;

correction = 0;

t_strt = 0.1;

t_end = 10000;

fitCount_tot = 5;

pause_num_array = [2, 3, 4, 5];

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

DWT_folder = [output_foldername slash 'selected_DWT_pdf_Corrected' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_DWT = [DWT_folder slash 'DWT_pdf_Corrected' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_DWT_dat = [DWT_folder slash 'DWT_Corrected' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

DWT_hist_dat = load([str_DWT, '.txt'], '-ascii');

bins_exp = transpose(DWT_hist_dat(:,1));

bar_pos_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,2));

DWT_hist_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,3));

Err_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,4));

DWT_array = transpose(load([str_DWT_dat, '.txt'], '-ascii'));

DWT_array = DWT_array(DWT_array > t_strt & DWT_array < t_end);

[bar_pos, bins]=make_bins(Ti,bins_per_decade,number_of_decades,t0,correction);

[DWT_hist]=DwellTimeHist_v3(DWT_array, t0, bins);

dat_length = length(DWT_array);

binwidth = bins(2:length(bins)) - bins(1:(length(bins)-1));

Weights = dat_length * DWT_hist .* binwidth;

bool_vec = Weights > 0;

Weights = Weights(bool_vec);

bincenter = bar_pos(bool_vec);

for n=1:length(pause_num_array)
    
    pause_num = pause_num_array(n)
    
    gamma_nums = [Nw, ones([1,pause_num - 1])];

    ub = [-1.27   3.0*ones([1,pause_num]) 2.0*ones([1,pause_num])    0.6    0.36];

    lb = [-1.28    0.4*ones([1,2*pause_num])    0.2    0.32];
    
    params_0 = 0.5*(ub+lb);

    options = saoptimset('TolFun',10^-10, 'MaxFunEvals', 10^5);

    MLfun= @(params) ML_func_CG(params, B, Weights, bincenter, gamma_nums, Nw, t_strt, t_end, Dobcktrck);
    
    for j=1:fitCount_tot
        
        disp(j)

        [fit,fval] = simulannealbnd(MLfun, params_0, lb, ub, options);
        
        BIC = 2*(fval + (pause_num*log(dat_length)));
        
        [q_array, k_array, sigma_ln, t_bar] = MLparams2modelparams(fit, pause_num);
        
        fit_dat = zeros([length(fit), 3]);

        fit_dat(:,1) = transpose([k_array, q_array, t_bar, sigma_ln]);

        fit_dat(:,2) = transpose(fit);

        fit_dat(1,3) = fval;
        
        fit_dat(2,3) = BIC;
        
        T_array=10.^(-1:0.01:4);

    	y = Calculate_Model_v2(T_array, q_array, k_array, gamma_nums, sigma_ln, t_bar, Nw, B, t_strt, t_end, Dobcktrck);

        figure()

        ax = axes();

        PL = loglog(T_array, y, 'g', 'Linewidth', 2);

        set(ax,'XScale','log','YScale','log')

        hold on

        errorbar(bar_pos_exp ,DWT_hist_exp, Err_exp,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

        hold off
    
        ylim([10^-8, 5])

        str_title = ['Pause number = '  num2str(pause_num) ' ' 'Ts = '  num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

        title(str_title)

        xlabel('Dwell Time (s)')

        ylabel('PDF')
        
        str_fit = [Fit_folder slash 'Fit_Result' '_' 'Pause_num=' num2str(pause_num) '_' 'fit_num=' num2str(j) '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

        save([str_fit, '.txt'], 'fit_dat', '-ascii')

        saveas(PL, [str_fit, '.jpg'], 'jpg')

        saveas(PL, [str_fit, '.fig'], 'fig')
        
        close
               
    end
    
end
