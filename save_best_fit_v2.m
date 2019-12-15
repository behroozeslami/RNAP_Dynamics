cam_freq = 25;

t0 = 1/cam_freq;

f0 = 12;

f = 7.5;

B = exp(-f/f0);

Trace_time = 10^4;

Trace_num = 40;

W = 12;

Ts = (2*W+1);

T = Ts;

k = 0;

Nw = 10;

Pause_threshold = 100.64;

Pause_threshold_V = 100.64;

Vmax = 60;

Vmin = -40;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

Dobcktrck = 0;

pause_num = 4;

t_strt = 0.1;

t_end = 1200;

mainfoldername = 'Export_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

Best_Fit_folder = [output_foldername slash 'Best_Fit' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

if ~(exist(Best_Fit_folder,'dir')==7)
    
    mkdir(Best_Fit_folder)
    
end

Fit_folder = [output_foldername slash 'Fit_Results_BS' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_best_fit = [Fit_folder slash 'Best_Pause_Fit_ConfInt' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

Conf_int_dat = load([str_best_fit, '.txt'], '-ascii');

len = pause_num;

pause_k_array = Conf_int_dat(1,2:(len + 1));

Low_Err_pause_k_array = Conf_int_dat(2,2:(len + 1));

Up_Err_pause_k_array = Conf_int_dat(3,2:(len + 1));

q_array = Conf_int_dat(1,(len + 2):(2*len + 1));

Low_Err_q_array = Conf_int_dat(2,(len + 2):(2*len + 1));

Up_Err_q_array = Conf_int_dat(3,(len + 2):(2*len + 1));

t_bar = Conf_int_dat(1,(2*len + 2));

sigma_ln = Conf_int_dat(1,(2*len + 3));

Fit_kel_folder = [output_foldername slash 'kel_Fit_Results' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_kelFit_dat = [Fit_kel_folder slash 'kel_Fit_Results' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

kelFit_dat = load([str_kelFit_dat, '.txt'], '-ascii');

kel = kelFit_dat(1,4);

Err_kel = kelFit_dat(2,4);

k_array = [kel, pause_k_array];

Low_Err_k_array = [Err_kel, Low_Err_pause_k_array];

Up_Err_k_array = [Err_kel, Up_Err_pause_k_array];

pause_num = length(q_array);

gamma_nums = [Nw, ones([1,pause_num - 1])];

noise_V_pdf_folder = [output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise_V_dat = [noise_V_pdf_folder slash 'NoiseVelocity' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

noise_velocities_sample = load([str_noise_V_dat '.txt'],'-ascii');

V_pdf_folder = [output_foldername slash 'Velocity_pdf_NoLongPause' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_V_pdf = [V_pdf_folder slash 'PauseFreeVelocity_PDF' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

V_hist_dat = load([str_V_pdf, '.txt'], '-ascii');

velocity_pdf = V_hist_dat(:,2)';

binwidth = V_hist_dat(1,4);

DWT_folder = [output_foldername slash 'selected_DWT_pdf_Corrected' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_DWT = [DWT_folder slash 'DWT_pdf_Corrected' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

DWT_hist_dat = load([str_DWT, '.txt'], '-ascii');

bins_exp = transpose(DWT_hist_dat(:,1));

bar_pos_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,2));

DWT_hist_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,3));

Err_exp = transpose(DWT_hist_dat(1:length(bins_exp)-1,4));

[centers, simul_velocity_pdf] = Simulate_velocity_pdf_nolongpause_v2(k_array,q_array,Dobcktrck,t0,f,Trace_time,Trace_num,W,k,Nw,T,Pause_threshold_V,noise_velocities_sample,binwidth,Vmin,Vmax);

T_array=10.^(-1:0.01:3.5);

y = Calculate_Model_v2(T_array, q_array, k_array, gamma_nums, sigma_ln, t_bar, Nw, B, t_strt, t_end, Dobcktrck);

[Simul_DWT_hist] = Simulate_dwellTime_pdf(k_array,q_array,Dobcktrck,t0,B,Trace_time,Trace_num,0,0,Nw,0,bins_exp);

Simul_DWT_hist(bins_exp < (Ts*t0)) = 0;

k_array_dat = [k_array; Low_Err_k_array; Up_Err_k_array];

q_array_dat = [q_array; Low_Err_q_array; Up_Err_q_array];

figure()

ax = axes();

PL1 = plot(centers, velocity_pdf, 'rO', 'MarkerFacecolor', 'r');

hold on

plot(centers, simul_velocity_pdf, 'g', 'LineWidth', 2)

hold off

set(ax,'XScale','linear','YScale','log')
    
xlabel('Pause Free Velocity (bp/s)')

ylabel('PDF')

legend('Data', 'ML Fit (with simulation)', 'Location', 'South')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'T = ' num2str(T) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

figure()

ax = axes();

PL1_2 = plot(centers, velocity_pdf, 'rO', 'MarkerFacecolor', 'r');

hold on

plot(centers, simul_velocity_pdf, 'g', 'LineWidth', 2)

hold off

set(ax,'XScale','linear','YScale','linear')
    
xlabel('Pause Free Velocity (bp/s)')

ylabel('PDF')

legend('Data', 'ML Fit (with simulation)', 'Location', 'South')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'T = ' num2str(T) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

figure()

ax = axes();

PL2 = loglog(T_array, y, 'g', 'Linewidth', 2);

set(ax,'XScale','log','YScale','log')

hold on

errorbar(bar_pos_exp ,DWT_hist_exp, Err_exp,'rO','MarkerSize',6, 'MarkerFaceColor', 'r');

loglog(bar_pos_exp, Simul_DWT_hist, 'k--', 'LineWidth', 1.5)
    
hold off

ylim([10^-8, 5])

legend('ML Fit', 'Data', 'Simulation (ideal model)')

str_title = ['Ts =' ' ' num2str(2*W+1) ' ' 'Nw = ' num2str(Nw) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

xlabel('Dwell Time (s)')

ylabel('PDF')

Col = {'r','g','b','m','k','c','y'};

figure()

leg_k ={'k_{el}'};

PL3 = semilogy(1, k_array(1), 'O', 'Color', Col{1}, 'MarkerSize',6, 'MarkerFaceColor', Col{1});

hold on

for j=1:length(k_array)-1
    
    semilogy(j+1, k_array(j+1), 'O', 'Color', Col{j+1}, 'MarkerSize',6, 'MarkerFaceColor', Col{j+1})
    
    leg_k{j+1} = ['k_{' num2str(j) '}'];
    
end

legend(leg_k)

for j=1:length(k_array)
    
    errorbar(j, k_array(j), Low_Err_k_array(j), Up_Err_k_array(j), '.', 'Color', Col{j})
    
end

set(gca,'XTick',1:length(k_array))

set(gca,'XTickLabel',{})

xlim([0.5,length(k_array)+0.5])

ylabel('rates (bp/s)')

hold off

figure()

for j=1:length(q_array)
    
    PL4 = semilogy(j, q_array(j), 'O', 'Color', Col{j}, 'MarkerSize',6, 'MarkerFaceColor', Col{j});
    
    leg_q{j} = ['P_{' num2str(j) '}'];
    
    hold on
    
end

legend(leg_q)

for j=1:length(q_array)
    
    errorbar(j, q_array(j), Low_Err_q_array(j), Up_Err_q_array(j), '.', 'Color', Col{j})
    
end

set(gca,'XTick',1:length(q_array))

set(gca,'XTickLabel',{})

xlim([0.5,length(q_array)+0.5])

ylabel('pause probabilities')

hold off

str_Fit_V_pdf_log = [Best_Fit_folder slash 'BestFit_PauseFreeVelocity_PDF_log' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

saveas(PL1, [str_Fit_V_pdf_log , '.jpg'], 'jpg')

saveas(PL1, [str_Fit_V_pdf_log , '.fig'], 'fig')

str_Fit_V_pdf = [Best_Fit_folder slash 'BestFit_PauseFreeVelocity_PDF' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

saveas(PL1_2, [str_Fit_V_pdf , '.jpg'], 'jpg')

saveas(PL1_2, [str_Fit_V_pdf , '.fig'], 'fig')

str_best_pause_fit = [Best_Fit_folder slash 'BestFit_DwellTime_PDF' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

saveas(PL2, [str_best_pause_fit, '.jpg'], 'jpg')

saveas(PL2, [str_best_pause_fit, '.fig'], 'fig')

str_k_array = [Best_Fit_folder slash 'BestFit_k_array' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'T=' num2str(T) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

saveas(PL3, [str_k_array, '.jpg'], 'jpg')

saveas(PL3, [str_k_array, '.fig'], 'fig')

save([str_k_array, '.txt'], 'k_array_dat', '-ascii')

str_q_array = [Best_Fit_folder slash 'BestFit_q_array' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'T=' num2str(T) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

saveas(PL4, [str_q_array, '.jpg'], 'jpg')

saveas(PL4, [str_q_array, '.fig'], 'fig')

save([str_q_array, '.txt'], 'q_array_dat', '-ascii')