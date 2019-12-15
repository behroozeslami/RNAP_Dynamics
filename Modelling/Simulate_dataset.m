f = 7.5;

force = abs(f);

f0 = 12;

delta = 0.7;

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0=1/cam_freq;

mainfoldername = 'NEW_AF_7.5pN_1mM_ITP_100uM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end


if ~(exist(output_foldername,'dir')==7)
    
    mkdir(output_foldername)
    
end

directory = ['../..' slash 'Richard_Data' slash mainfoldername slash subfolder];

if strcmp('',subfolder)
    directory = ['../..' slash 'Richard_Data' slash mainfoldername];
end

WLC=xlsread(['../..' slash 'Richard_Data' slash 'WLC.xlsx']);

force_array = WLC(:,1);

bp_step_array = WLC(:,2);

bp_step_Modelfit = mean(bp_step_array(abs(force_array - force) == min(abs(force_array - force))));

info = dir(directory);

filenames = {info.name};

filenames = filenames(3:length(filenames));

trace_num = length(filenames);

str_Param = 'Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx';

[Param_dat, txt] = xlsread(str_Param);

Params = Param_dat(1,:);

% Parameters for GreB OF

force_indices = [6];

kb = Params(force_indices(1));

%Params(force_indices(1)) = kb*exp((7.5-f)*(1-delta)/f0);

Params(1) = 19.6;

Params(2) = 0.194;

Params(6) = 0.732;

% Noise data

str_noise = '../OF4_9pN_1mM_NTP_GreB_M&H/Noise_Velocity_pdf_Ts=49_T=49_Nw=4_Pause_threshold_V=100.64s/NoiseSample_small_Ts=49_l0_range=0.13-0.46_Seidel_range=0-0.5_Nw=4_Pause_threshold_V=100.64s.mat';

noise_dat = load(str_noise, '-mat');
    
noise = transpose(noise_dat.noise_dat_small);

output_foldername_simul = ['../..' slash 'Richard_Data' slash 'Sim_AF_7.5pN_1mM_ITP_100uM_NTP'];

if ~(exist(output_foldername_simul,'dir')==7)
    
    mkdir(output_foldername_simul)
    
end

for n=1:trace_num
    
    filename = filenames(n);
    
    str_load = [directory slash filename{1}];

    dat = load(str_load , '-ascii');

    t_array = dat(:,1);
    
    Trace_time = max(t_array)+t0;
    
    bp_step_Exp = dat(1,10);
    
    [t_array, x_array] = Gillespie_MicroModel_1(Params, t0, Trace_time);
    
    points = length(x_array);

    x_noise = gen_Exp_noise_v3(noise, points);

    x_array = x_array + x_noise;

    simul_dat = dat;
    
    simul_dat(:,8) = x_array(1:length(dat))';

    save([output_foldername_simul slash filename{1}], 'simul_dat', '-ascii')
    
end


