f0 = 12;

f = 7.5;

B = exp(-f/f0);

cam_freq = 25;

t0 = 1/cam_freq; 

Trace_time = 10^4;

Trace_num = 100;

W = 12;

Ts = (2*W+1);

T = Ts;

k = 0;

Nw = 10;

Pause_threshold = 100.64;

Pause_threshold_V = 100.64;

Vmax = 60;

Vmin = -40;

Vmax_fit = 35;

Vmin_fit = -20;

Fit_Count_tot = 100;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

Dobcktrck = 0;

t_strt = 0.1;

t_end = 1200;

mainfoldername = 'Export_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

Fit_kel_folder = [output_foldername slash 'kel_Fit_Results' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

if ~(exist(Fit_kel_folder,'dir')==7)
    
    mkdir(Fit_kel_folder)
    
end

kel_array = 24:0.05:26;

Fit_folder = [output_foldername slash 'Fit_Results_BS' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

str_best_fit = [Fit_folder slash 'Best_Pause_Fit' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'BT=' num2str(Dobcktrck) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Time_range=' num2str(t_strt) '-' num2str(t_end) '_' 'Pause_threshold=' num2str(Pause_threshold) 's'];

best_model_param_dat = load([str_best_fit, '.txt'], '-ascii');

Size = size(best_model_param_dat);

len = Size(2);

pause_k_array = best_model_param_dat(1,2:len);

q_array = best_model_param_dat(3,1:(len-1));

%{
pause_k_array = [0.5474    0.0418    0.0021];

q_array = [0.0664    0.0078    0.0014];
%}

noise_V_pdf_folder = [output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise_V_dat = [noise_V_pdf_folder slash 'NoiseVelocity' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

noise_velocities_sample = load([str_noise_V_dat '.txt'],'-ascii');

V_pdf_folder = [output_foldername slash 'Velocity_pdf_NoLongPause' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_V_pdf = [V_pdf_folder slash 'PauseFreeVelocity_PDF' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

V_hist_dat = load([str_V_pdf, '.txt'], '-ascii');

centers = V_hist_dat(:,1)';

velocity_pdf = V_hist_dat(:,2)';

Err = V_hist_dat(:,3)';

binwidth = V_hist_dat(1,4);

dat_length = V_hist_dat(2,4);

Weights = (1/Ts)*dat_length * velocity_pdf * binwidth;

Weights(centers > Vmax_fit) = 0;

Weights(centers < Vmin_fit) = 0;

sum_ML_func_array = zeros([1, length(kel_array)]);

SQsum_ML_func_array = zeros([1, length(kel_array)]);

ML_func_array_now = zeros([1, length(kel_array)]);

sum_best_kel = 0;

for Fit_Count=1:Fit_Count_tot
    
    for j=1:length(kel_array)
    
        k_array = [kel_array(j), pause_k_array];
    
        [centers, simul_velocity_pdf] = Simulate_velocity_pdf_nolongpause_v2(k_array,q_array,Dobcktrck,t0,B,Trace_time,Trace_num,W,k,Nw,T,Pause_threshold_V,noise_velocities_sample,binwidth,Vmin,Vmax);
        
        log_simul_velocity_pdf = log(simul_velocity_pdf);
        
        log_simul_velocity_pdf(centers > Vmax_fit) = 0;

        log_simul_velocity_pdf(centers < Vmin_fit) = 0;
        
        ML_func_array_now(j) = -dot(Weights, log_simul_velocity_pdf);
    
        sum_ML_func_array(j) = sum_ML_func_array(j) + ML_func_array_now(j) ;
        
        SQsum_ML_func_array(j) = SQsum_ML_func_array(j) + ML_func_array_now(j)^2;
        
        disp(Fit_Count)
    
        disp(kel_array(j))
    
    end
    
    [Min_ML_func_now, Index] = min(ML_func_array_now);
    
    best_kel_now = kel_array(Index);
    
    sum_best_kel = sum_best_kel + best_kel_now;
    
end

ML_func_array = sum_ML_func_array/Fit_Count;

Err_ML_func_array = sqrt((SQsum_ML_func_array/Fit_Count) - ML_func_array.^2)/sqrt(Fit_Count-1);

best_kel = sum_best_kel/Fit_Count;

[Min, Index_0] = min(abs(kel_array - best_kel));

Min_ML_func = ML_func_array(Index_0);

y = transpose(ML_func_array - Min_ML_func);

x = transpose(kel_array - best_kel);

parabola = fittype('deltay + 0.5*(x-deltax)^2/sigma^2');

[parabola_fit,gof,fitinfo] = fit(x,y,parabola,'Startpoint',[0,0,0.1]);

x = x - parabola_fit.deltax;

y = y - parabola_fit.deltay;

best_kel = best_kel + parabola_fit.deltax;

Err_kel = parabola_fit.sigma;

Fit_dat = [kel_array', y, Err_ML_func_array', [best_kel, Err_kel, zeros([1, length(kel_array)-2])]'];

k_array = [best_kel, pause_k_array];

[centers, simul_velocity_pdf] = Simulate_velocity_pdf_nolongpause_v2(k_array,q_array,Dobcktrck,t0,B,Trace_time,Trace_num,W,k,Nw,T,Pause_threshold_V,noise_velocities_sample,binwidth,Vmin,Vmax);

figure()

ax = axes();

PL = plot(centers, velocity_pdf, 'rO', 'MarkerFacecolor', 'r');

hold on

plot(centers, simul_velocity_pdf, 'g', 'LineWidth', 2)

line([Vmax_fit, Vmax_fit], ylim)

line([Vmin_fit, Vmin_fit], ylim)

hold off

set(ax,'XScale','linear','YScale','linear')
    
xlabel('Pause Free Velocity (bp/s)')

ylabel('PDF')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'T = ' num2str(T) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

figure()

PL2 = errorbar(kel_array, y', Err_ML_func_array, 'rO', 'MarkerFacecolor', 'r');

hold on

plot(kel_array,0.5*x'.^2/Err_kel^2,'k','LineWidth', 2)

hold off

xlabel('Elongation rate (bp/s)')

ylabel('ML Func')

str_Fit_V_pdf = [Fit_kel_folder slash 'PauseFreeVelocity_PDF' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_Fit_V = [Fit_kel_folder slash 'PauseFreeVelocity' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_Fit_dat = [Fit_kel_folder slash 'kel_Fit_Results' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

saveas(PL, [str_Fit_V_pdf , '.jpg'], 'jpg')

saveas(PL, [str_Fit_V_pdf , '.fig'], 'fig')

saveas(PL2, [str_Fit_dat , '.jpg'], 'jpg')

saveas(PL2, [str_Fit_dat , '.fig'], 'fig')

save([str_Fit_dat, '.txt'], 'Fit_dat', '-ascii')

