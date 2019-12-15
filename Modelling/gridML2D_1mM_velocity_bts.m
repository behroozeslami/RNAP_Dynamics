function [] = gridML2D_1mM_velocity_bts(kel, P)

disp('Started at')

disp(clock)

if ischar(kel)
    
    kel = str2double((kel));
    
end

if ischar(P)
    
    P = str2double((P));
    
end

Trace_time = 3*10^3;

V_num = 10^5;

fit_index_1 = 1;

fit_index_2 = 2;

Vmin_fit = -20;

Vmax_fit = 30;

Vmax = 60;

Vmin = -40;

binwidth = 1;

W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 10;

Pause_threshold_V = 100.64;

cam_freq = 25;

t0 = 1/cam_freq;

Model_ID = 1;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'Select_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

if ~(exist(output_foldername,'dir')==7)
    
    mkdir(output_foldername)
    
end

str_V = [output_foldername '/V_Bst' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T) '_' 'limits=' num2str(Vmin_fit) '-' num2str(Vmax_fit) '.mat'];

D = load(str_V, '-mat');

V_Mat = D.V_Mat; 

noise_V_pdf_folder = ['..' slash output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise = [noise_V_pdf_folder slash 'NoiseSample_small' '_' 'Ts='  num2str(2*W+1) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '.mat'];

noise_dat = load(str_noise, '-mat');
    
noise = transpose(noise_dat.noise_dat_small);

noise_type = 'Exp';

[Param_dat, txt] = xlsread('Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Params = Param_dat(1,:);

Param_name_1 = txt(1,fit_index_1+1);

Param_name_2 = txt(1,fit_index_2+1);

Fit_folder = [output_foldername slash 'V_Fit_Results' '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T*t0) 's'];

if ~(exist(Fit_folder,'dir')==7)
    
    mkdir(Fit_folder)
    
end

str_Fit_dat = [Fit_folder slash 'V_Fit_Results' '_' Param_name_1{1} '=' num2str(kel) '_', Param_name_2{1} '=' num2str(P) '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T*t0) 's' '_' 'limits=' num2str(Vmin_fit) '-' num2str(Vmax_fit)];

str_simul_dat = [Fit_folder slash 'Simul_Vpdf' '_' Param_name_1{1} '=' num2str(kel) '_', Param_name_2{1} '=' num2str(P) '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T*t0) 's'];

Params(fit_index_1) = kel; 

Params(fit_index_2) = P;

disp('Parameters')

disp(Params)

disp('Simulation Started!')

[Fit_dat, simul_velocity_pdf] = Calculate_Velocity_LogLikelihood_bootstrapping(Params, V_Mat, Vmin_fit, Vmax_fit, Model_ID, t0, Trace_time, V_num, W, T, Nw, Pause_threshold_V, noise, noise_type, Vmin, Vmax, binwidth);

data = [kel, P, Fit_dat];

Simul_data = [kel, P, simul_velocity_pdf];

save([str_Fit_dat, '.txt'], 'data', '-ascii')

save([str_simul_dat, '.txt'], 'Simul_data', '-ascii')

disp('Finished at')

disp(clock)

exit


