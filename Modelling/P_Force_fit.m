P_array = zeros([1,9]);

F_array = [-12.5, -10, -9, -7.5, -5, 5, 7.5, 10 12.5];

LE_array = zeros([1, 9]);

UE_array = zeros([1, 9]);

EM_array = zeros([1, 9]);

slash = '/';

% 12AF

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Nw = 4;

t_strt = 0.05;

t_end = 10000;

output_foldername = 'Export_AF_12.5pN_1mM_NTP';

save_folder = ['Analysis_12AF_delta=0.7_SingleTraces_kel-P' slash '12AF_S'];

data_name = '12AF_S';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(9) = Fit_dat(1,2);

LE_array(9) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(9) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(9) = Fit_dat(4,2);

% 10 AF

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Nw = 4;

t_strt = 0.05;

t_end = 10000;

output_foldername = 'EXPORT_AF_10pN_1mM_NTP';

save_folder = ['Analysis_10AF_delta=0.7_SingleTraces_kel-P' slash '10AF_HV2'];

data_name = '10AF_HV2';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(8) = Fit_dat(1,2);

LE_array(8) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(8) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(8) = Fit_dat(4,2);

% 7.5AF

cam_freq = 25; 

t0 = 1/cam_freq;

W = 12;

Nw = 4;

t_strt = 0.05;

t_end = 10000;

output_foldername = 'Export_AF_7.5pN_1mM_NTP';

save_folder = ['Analysis_1mM_SingleTraces_kel-P' slash '1mM_MV_2'];

data_name = '1mM_MV_2';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(7) = Fit_dat(1,2);

LE_array(7) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(7) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(7) = Fit_dat(4,2);

% 5 AF

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Nw = 4;

t_strt = 0.05;

t_end = 10000;

output_foldername = 'EXPORT_AF_5pN_1mM_NTP';

save_folder = ['Analysis_5AF_delta=0.7_SingleTraces_kel-P' slash '5AF_MV_1'];

data_name = '5AF_MV_1';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(6) = Fit_dat(1,2);

LE_array(6) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(6) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(6) = Fit_dat(4,2);

% 5 OF

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Nw = 4;

t_strt = 0.05;

t_end = 10000;

output_foldername = 'EXPORT_OF_5pN_1mM_NTP';

save_folder = ['Analysis_5OF_delta=0.7_SingleTraces_kel-P' slash '5OF_MV_1'];

data_name = '5OF_MV_1';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(5) = Fit_dat(1,2);

LE_array(5) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(5) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(5) = Fit_dat(4,2);

% 7.5 OF

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Nw = 4;

t_strt = 0.05;

t_end = 10000;

output_foldername = 'EXPORT_OF_7.5pN_1mM_NTP_NoB2';

save_folder = ['Analysis_7OF_delta=0.7_SingleTraces_kel-P' slash '7.5OF_MV1'];

data_name = '7.5OF_MV1';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(4) = Fit_dat(1,2);

LE_array(4) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(4) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(4) = Fit_dat(4,2);

% 9 OF

cam_freq = 50;

t0 = 1/cam_freq;

W = 24;

Nw = 4;

t_strt = 0.05;

t_end = 10000;

output_foldername = 'OF4_9pN_1mM_NTP_Control';

save_folder = ['Analysis_9OF_delta=0.7_SingleTraces_kel-P' slash '9OF_MV_1'];

data_name = '9OF_MV_1';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(3) = Fit_dat(1,2);

LE_array(3) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(3) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(3) = Fit_dat(4,2);

% 10 OF

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Nw = 4;

t_strt = 0.05;

t_end = 10000;

output_foldername = 'EXPORT_OF_10pN_1mM_NTP';

save_folder = ['Analysis_10OF_delta=0.7_SingleTraces_kel-P' slash '10OF_MV_2'];

data_name = '10OF_MV_2';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(2) = Fit_dat(1,2);

LE_array(2) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(2) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(2) = Fit_dat(4,2);

% 12 OF

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Nw = 4;

t_strt = 0.05;

t_end = 10000;

output_foldername = 'EXPORT_OF_12.5pN_1mM_NTP';

save_folder = ['Analysis_12OF_delta=0.7_SingleTraces_kel-P' slash '12OF_S_0.7'];

data_name = '12OF_S_0.7';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(1) = Fit_dat(1,2);

LE_array(1) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(1) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(1) = Fit_dat(4,2);


x = F_array';

y = P_array';

PM_ref = 0.95; % inferred from the NTP fit

K_ref = 14; % inferred from the NTP fit

func =@(delta,x) 1-(1+(1/PM_ref-1)*exp(-(x-7.5)*delta/12)).^(-1)./(1+0.001*K_ref*exp((7.5-x)*(1-delta)/12));

P_F = fittype(func);

PF_fit = fit(x,y,P_F, 'Lower',[0], 'Upper',[1], 'Start', [0.5]);

PL = plot(PF_fit,'b-');

hold on 

errorbar(F_array, P_array, LE_array, UE_array, 'rO','MarkerSize',6, 'MarkerFaceColor', 'r')

errorbar(F_array, P_array, EM_array, 'k.')

xlabel('force (pN)')

ylabel('P_{tot}')

hold off


folder_name = 'Force_dependency';

if ~(exist(folder_name,'dir')==7)
    
    mkdir(folder_name)
    
end

str_save = [folder_name '/' 'P-F'];

saveas(PL,[str_save,'.fig'], 'fig')

saveas(PL,[str_save,'.jpg'], 'jpg')


P_F_dat = zeros([5, length(F_array)]);

P_F_dat(1,:) = F_array;

P_F_dat(2,:) = P_array;

P_F_dat(3,:) = LE_array;

P_F_dat(4,:) = UE_array;

P_F_dat(5,:)= EM_array;

row_header(1:5,1) = {'Force (pN)', 'P_tot', 'LE', 'UE', 'std_of_mean'};

xlswrite([str_save '.xlsx'],P_F_dat,'Sheet1','B1');

xlswrite([str_save '.xlsx'],row_header,'Sheet1','A1'); 

F_array_large = (-13:0.1:13)';

FIT = func(PF_fit.delta, F_array_large);

Theory = [F_array_large, FIT];

col_names = {'Force (pN)', 'P_tot (theory)'};

xlswrite([str_save '.xlsx'],Theory,'Sheet2','A2');

xlswrite([str_save '.xlsx'],col_names,'Sheet2','A1'); 





