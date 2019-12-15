function [] = gridML2D_12OF_BCKTR_W_delta(Wall, delta)

f = -12.5;

kel = 18.5;

P = 0.29;

disp('Started at')

disp(clock)

if ischar(Wall)
    
    Wall = str2double((Wall));
    
end

if ischar(delta)
    
    delta = str2double((delta));
    
end

Trace_time = 2100;

DWT_num = 10^5;

t_strt = 0.05;

t_end = 10000;

W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 4;

Pause_threshold_V = 100.64;

cam_freq = 25;

t0 = 1/cam_freq;

Model_ID = 2;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'Select_OF_12.5pN_1mM_NTP';

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

str_DWT_dat = [DWT_folder slash 'DWT_Combined' '_' 'Ts='  num2str(2*W+1) '_' 'Nw=' num2str(Nw) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

DWT_Mat = transpose(load([str_DWT_dat, '.txt'], '-ascii'));

dat_length = length(DWT_Mat);

noise_V_pdf_folder = ['..' slash output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise = [noise_V_pdf_folder slash 'NoiseSample_small' '_' 'Ts='  num2str(2*W+1) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '.mat'];

noise_dat = load(str_noise, '-mat');
    
noise = transpose(noise_dat.noise_dat_small);

noise_type = 'Exp';

[Param_dat, txt] = xlsread('Parameters/Backtrack_Model_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Params = Param_dat(1,:);

Param_name_1 = {'W'};

Param_name_2 = {'delta'};

Fit_folder = [output_foldername slash 'Fit_Results' '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

if ~(exist(Fit_folder,'dir')==7)
    
    mkdir(Fit_folder)
    
end

str_Fit_dat = [Fit_folder slash 'Fit_Results' '_' Param_name_1{1} '=' num2str(Wall) '_', Param_name_2{1} '=' num2str(delta) '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end)];

str_simul_dat = [Fit_folder slash 'Simul_DWTpdf' '_' Param_name_1{1} '=' num2str(Wall) '_', Param_name_2{1} '=' num2str(delta) '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

str_proc = [Fit_folder slash 'Processivity' '_' Param_name_1{1} '=' num2str(Wall) '_', Param_name_2{1} '=' num2str(delta)];

Params(1) = kel;

Params(2) = P;

Params = [Params, delta, f, Wall];

disp('Parameters')

disp(Params)

disp('Simulation Started!')

[Fit_dat, Simul_DWT_array, Proc_data] = Calculate_LogLikelihood_SingleTrace(Params, dat_length, DWT_Mat, t_strt, t_end, Model_ID, t0, Trace_time, DWT_num, W, Nw, noise, noise_type);

Simul_DWT_hist = DwellTimeHist_v3(Simul_DWT_array, t0, bins_exp);

data = [Wall, delta, Fit_dat];

Simul_data = [Wall, delta, Simul_DWT_hist];

Proc = [Wall, delta, Proc_data];

save([str_Fit_dat, '.txt'], 'data', '-ascii')

save([str_simul_dat, '.txt'], 'Simul_data', '-ascii')

save([str_proc, '.txt'], 'Proc', '-ascii')

disp('Finished at')

disp(clock)

exit


