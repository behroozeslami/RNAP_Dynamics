W = 12;

k = 0;

Nw = 4;

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0 = 1/cam_freq;

t_strt = 0.05;

t_end = 10000;

mainfoldername = 'NEW_AF_7.5pN_1mM_ITP_100uM_NTP';

data_name = 'ITP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

DWT_folder = ['Modelling' slash output_foldername];

if ~(exist(DWT_folder,'dir')==7)
    
    mkdir(DWT_folder)
    
end

directory = ['..' slash 'Richard_Data' slash mainfoldername slash subfolder];

if strcmp('',subfolder)
    
    directory = ['..' slash 'Richard_Data' slash mainfoldername];
    
end

str_trace_index = [DWT_folder slash 'Trace_indices' '.mat'];

trace_indices = load(str_trace_index);

trace_indices = trace_indices.selected_trace_number;

str_DWT_length = [DWT_folder slash 'DWT_data_length' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

dat_length = load(str_DWT_length);

DWT_length = dat_length.dat_length;

str_DWT_dat = [DWT_folder slash 'DWT_data' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

DWTs = load(str_DWT_dat);

DWTs = DWTs.DWT_dat;

str_V = [DWT_folder slash 'EndtoEnd_V' '.mat'];

V = load(str_V);

V = V.V_array;

Max_DWT = arrayfun(@(row) max(DWTs(row, :)), 1:length(trace_indices));

opt = statset('MaxIter',1000);

GMM_V2 = fitgmdist(V',2,'Options', opt);

GMM_V3 = fitgmdist(V',3,'Options', opt);

GMM_V = GMM_V2;

GV_num = 2;

V_names ={'LV', 'HV'};

if GMM_V3.AIC < GMM_V2.AIC
    
    GMM_V = GMM_V3;
    
    GV_num = 3;
    
    V_names ={'LV', 'MV', 'HV'};
    
end

V_mu = GMM_V.mu;

[sorted_V_mu, I] = sort(V_mu);

ids_V = 1:GV_num;

ids_V = ids_V(I);

idx = cluster(GMM_V, V'); 

for j=1:GV_num
    
    logMax_DWT = log(Max_DWT(idx==ids_V(j)));
    
    selected_trace_indices_V = trace_indices(idx==ids_V(j));
    
    GMM_Max_DWT1 = fitgmdist(logMax_DWT',1,'Options', opt);
    
    GMM_Max_DWT2 = fitgmdist(logMax_DWT',2,'Options', opt);
    
    GMax_DWT_num = 1;
    
    GMM_Max_DWT = GMM_Max_DWT1;
    
    if GMM_Max_DWT2.AIC <GMM_Max_DWT1.AIC
    
        GMM_Max_DWT = GMM_Max_DWT2;
    
        GMax_DWT_num = 2;
    
    end
    %%%%%%%
    
    %GMM_Max_DWT = GMM_Max_DWT2;
    
    %GMax_DWT_num = 2;
     
    Max_DWT_mu = GMM_Max_DWT.mu;

    [sorted_Max_DWT_mu, I] = sort(Max_DWT_mu);

    ids_Max_DWT = 1:GMax_DWT_num;

    ids_Max_DWT = ids_Max_DWT(I);

    idx2 = cluster(GMM_Max_DWT, logMax_DWT'); 
    
    for n=1:GMax_DWT_num
        
        selected_trace_indices = selected_trace_indices_V(idx2==ids_Max_DWT(n))';
        
        V_name = V_names(j);
        
        class_name = [DWT_folder slash data_name, '_' , V_name{1}, '_', num2str(n) '.txt'];
        
        save(class_name, 'selected_trace_indices', '-ascii')
        
    end
    
end