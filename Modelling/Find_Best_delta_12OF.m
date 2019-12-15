f = -12.5;

f0 = 12;

fit_index_1 = 1;

fit_index_2 = 2;

param_num = 2;

t_strt = 0.05;

t_end = 10000;

new_t_strt = 0.05;

%new_t_end = 350;
new_t_end = t_end;

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 4;

Pause_threshold_V = 100.64; 

Trace_time = 3*10^3;

DWT_num = 10^8;

Model_ID = 1;

alpha = 0.136; % 1-sigma confidence interval

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'EXPORT_OF_12.5pN_1mM_NTP';

selected_trace_filename = '12OF_SM';

selected_trace_indices = load([mainfoldername '/' 'Selection_maxDWT_GMM3' '/' selected_trace_filename '.txt'])';

trace_indices = load([mainfoldername '/' 'Trace_indices.mat']);

trace_indices = trace_indices.selected_trace_number;

str_DWT_dat_single = [mainfoldername '/' 'DWT_data' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

D = load(str_DWT_dat_single, '-mat');

DWT_Mat = D.DWT_dat; 

DWT_Mat = DWT_Mat(ismember(trace_indices,selected_trace_indices),:);

DWT_array = reshape(DWT_Mat,[1, numel(DWT_Mat)]);

DWT_array = DWT_array(DWT_array>0);

DWT_array = DWT_array((DWT_array>=new_t_strt)&(DWT_array<=new_t_end));

noise_V_pdf_folder = ['..' '/' mainfoldername '/' 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise = [noise_V_pdf_folder '/' 'NoiseSample_small' '_' 'Ts='  num2str(2*W+1) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '.mat'];

noise_dat = load(str_noise, '-mat');
    
noise = transpose(noise_dat.noise_dat_small);

noise_type = 'Exp';

[Param_dat, txt] = xlsread('Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Params = Param_dat(1,:);

Param_name_1 = txt(1,fit_index_1+1);

Param_name_2 = txt(1,fit_index_2+1);

force_indices = [6];

kb = Params(force_indices(1));

delta_array = 0:0.1:1;

L = length(delta_array);

ML_score = zeros([1,L]);

for j=1:L
    
    delta = delta_array(j);
    
    disp(delta)
    
    subfolder_name = [selected_trace_filename '_' num2str(delta)];
    
    data_name = subfolder_name;
    
    fit_folder = ['Analysis_12OF_delta=' num2str(delta) '_SingleTraces_kel-P'];
    
    fit_subfolder = [mainfoldername '/' fit_folder '/' subfolder_name];
    
    str_FitResult = [fit_subfolder '/' 'Fit_outcome' '_' data_name '_' Param_name_1{1} '-' Param_name_2{1} '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];
    
    fit_data = xlsread(str_FitResult);
    
    kel = fit_data(1,1);
    
    P = fit_data(1,2);
    
    Params(fit_index_1) = kel; 

    Params(fit_index_2) = P;
    
    Params(force_indices(1)) = kb*exp((7.5-f)*(1-delta)/f0);
    
    ML_score(j) = Calculate_LogLikelihood_SingleTrace(Params, length(DWT_array), DWT_array, new_t_strt, new_t_end, Model_ID, t0, Trace_time, DWT_num, W, Nw, noise, noise_type);
    
end

ML_score = ML_score - min(ML_score);

delta_dat = [delta_array', ML_score'];

P_delta = exp(-ML_score);

P_delta = P_delta/sum(P_delta);

[B,M,L,U] = Conf_int(delta_array,P_delta,alpha);

delta_conf = [B,M,L,U]';

save_folder = [mainfoldername '/' 'Best_delta_12OF_SM'];

if ~(exist(save_folder,'dir')==7)
    
    mkdir(save_folder)
    
end

PL1 = figure();

bar(delta_array,P_delta,'r')

hold on

line([L,L],ylim,'LineWidth', 2, 'Color', 'b')

line([U,U],ylim,'LineWidth', 2, 'Color', 'b')

line([M,M],ylim,'LineWidth', 2, 'Color', 'g')

hold off

xlabel('delta')

ylabel('probability')

title('12OF')

xlim([0,1])

Pdelta_filename = [save_folder '/' 'Pdelta_12OF'];

saveas(PL1,[Pdelta_filename,'.fig'], 'fig')

saveas(PL1,[Pdelta_filename,'.jpg'], 'jpg')


col_header = {'delta'};

row_header(1:4,1) = {'Best', 'Mean', '1-sigma LB', '1-sigma UB'};

str_deltaFit = [save_folder '/'  'delta_Fit_outcome_12OF' '.xlsx'];

xlswrite(str_deltaFit,delta_conf,'Sheet1','B2');

xlswrite(str_deltaFit,col_header,'Sheet1','B1');

xlswrite(str_deltaFit,row_header,'Sheet1','A2'); 


figure()

PL2 = plot(delta_array, ML_score, 'rO--','MarkerSize',6, 'MarkerFaceColor', 'r');

xlabel('delta')

ylabel('LogLikelihood Score')

title('12OF')

delta_filename = [save_folder '/' 'LogLikelihood_vs_delta_12OF'];

saveas(PL2,[delta_filename,'.fig'], 'fig')

saveas(PL2,[delta_filename,'.jpg'], 'jpg')

save([delta_filename,'.txt'], 'delta_dat', '-ascii');


