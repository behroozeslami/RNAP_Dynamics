function [] = gridML2D_9OF_GreA_ke2_kp3_bts(ke2, kp3)

disp('Started at')

disp(clock)

f = -9;

f0 = 12;

delta = 0.7;

if ischar(ke2)
    
    ke2 = str2double((ke2));
    
end

if ischar(kp3)
    
    kp3 = str2double((kp3));
    
end

Trace_time = 3*10^3;

DWT_num = 10^5;

fit_index_1 = 4;

fit_index_2 = 7;

t_strt = 0.05;

t_end = 10000;

W = 24;

Ts = (2*W+1);

T = Ts;

Nw = 4;

Pause_threshold_V = 100.64;

cam_freq = 50;

t0 = 1/cam_freq;

Model_ID = 1;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

%mainfoldername = 'OF4_9pN_1mM_NTP_GreA_M';
mainfoldername = 'OF4_9pN_1mM_NTP_GreA_c';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

if ~(exist(output_foldername,'dir')==7)
    
    mkdir(output_foldername)
    
end

DWT_folder = ['..' slash output_foldername slash 'selected_DWT_pdf_Combined' '_' 'Ts=' num2str(2*W+1) '_' 'Nw=' num2str(Nw)];

str_DWT = [DWT_folder slash 'DWT_pdf_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

DWT_hist_dat = load([str_DWT, '.txt'], '-ascii');

bins_exp = transpose(DWT_hist_dat(:,1));

str_DWT_length = [output_foldername slash 'DWT_data_length_withBS' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

str_DWT_dat_single = [output_foldername slash 'DWT_data_withBS' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

DL = load(str_DWT_length, '-mat');

D = load(str_DWT_dat_single, '-mat');

DWT_Mat = D.DWT_dat; 

dat_length = DL.dat_length;

noise_V_pdf_folder = ['..' slash output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise = [noise_V_pdf_folder slash 'NoiseSample_small' '_' 'Ts='  num2str(2*W+1) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '.mat'];

noise_dat = load(str_noise, '-mat');
    
noise = transpose(noise_dat.noise_dat_small);

noise_type = 'Exp';

[Param_dat, txt] = xlsread('Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Params = Param_dat(1,:);

Param_name_1 = txt(1,fit_index_1+1);

Param_name_2 = txt(1,fit_index_2+1);

Fit_folder = [output_foldername slash 'Fit_Results' '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

if ~(exist(Fit_folder,'dir')==7)
    
    mkdir(Fit_folder)
    
end

str_Fit_dat = [Fit_folder slash 'Fit_Results' '_'  Param_name_1{1} '=' num2str(ke2) '_', Param_name_2{1} '=' num2str(kp3) '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

str_simul_dat = [Fit_folder slash 'Simul_DWTpdf' '_' Param_name_1{1} '=' num2str(ke2) '_', Param_name_2{1} '=' num2str(kp3) '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

Params(fit_index_1) = ke2; 

Params(fit_index_2) = kp3;

Params(1) = 21; % Control

Params(2) = 0.1; % Control

force_indices = [6];

kb = Params(force_indices(1));

Params(force_indices(1)) = kb*exp((7.5-f)*(1-delta)/f0);

disp('Parameters')

disp(Params)

disp('Simulation Started!')

[Fit_dat, Simul_DWT_array] = Calculate_LogLikelihood_SingleTrace(Params, dat_length, DWT_Mat, t_strt, t_end, Model_ID, t0, Trace_time, DWT_num, W, Nw, noise, noise_type);

Simul_DWT_hist = DwellTimeHist_v3(Simul_DWT_array, t0, bins_exp);

data = [ke2, kp3, Fit_dat];

Simul_data = [ke2, kp3, Simul_DWT_hist];

save([str_Fit_dat, '.txt'], 'data', '-ascii')

save([str_simul_dat, '.txt'], 'Simul_data', '-ascii')

disp('Finished at')

disp(clock)

exit


