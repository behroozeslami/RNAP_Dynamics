W = 12;

Ts = (2*W+1);

T = Ts;

k = 0;

Nw = 10;

%Pause_threshold_V = 100.64;

Pause_threshold_V = 5;

nboot=100;

cam_freq = 25;

t0 = 1/cam_freq;

binwidth = 1;

Vmax = 60;

Vmin = -40;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'heterogeneity';

subfolder = 'Fast';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

V_pdf_folder = [output_foldername slash 'Velocity_pdf_NoLongPause_v2' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

if ~(exist(V_pdf_folder,'dir')==7)
    
    mkdir(V_pdf_folder)
    
end

directory = ['..' slash 'Richard_Data' slash mainfoldername slash subfolder];

if strcmp('',subfolder)
    directory = ['..' slash 'Richard_Data' slash mainfoldername];
end

file_str_txt = [output_foldername slash 'Trace_data' '.txt'];

dat = load(file_str_txt, '-ascii');

trace_number = transpose(dat(:,1));

seidel_array = transpose(dat(:,2));

l0_array = transpose(dat(:,3));

l0_WLC = dat(1,4);

select_vec = (l0_array >= l0_min) & (l0_array <= l0_max) & (seidel_array >= Seidel_min) & (seidel_array <= Seidel_max);

selected_trace_number = floor(trace_number(select_vec));

selected_l0_array = l0_array(select_vec);

Blacklist = load([output_foldername '/' 'Blacklist.txt'], '-ascii');

High_V_Blacklist = load([output_foldername '/' 'Blacklist_High_Velocity.txt'], '-ascii');

Low_V_Blacklist = load([output_foldername '/' 'Blacklist_Low_Velocity.txt'], '-ascii');

Blacklist = [Blacklist; High_V_Blacklist; Low_V_Blacklist];

index_2_keep = Exclude_Blacklist(selected_trace_number, Blacklist);

selected_trace_number = selected_trace_number(index_2_keep);

selected_l0_array = selected_l0_array(index_2_keep);

str_V_pdf = [V_pdf_folder slash 'PauseFreeVelocity_PDF' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_V_dat = [V_pdf_folder slash 'PauseFreeVelocity' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

Dists = [];

for i=1:length(selected_trace_number)
    
    filename = ['RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = selected_l0_array(i);
    
    factor = l0/l0_WLC;

    x_array = factor*trace_dat(:,8);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);
    
    len = length(xs_array);
    
    Dist_array = transpose(xs_array(T+1:len) - xs_array(1:len-T));
    
    [Tc_array] = find_crossing_time(t_array, xs_array,Nw);
    
    [DWT_array]=find_Dwell_time(Tc_array);
    
    DWT_indices = 1:length(DWT_array);
    
    pause_start_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V));
    
    pause_end_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V)+1);
    
    pause_start_indices = round(pause_start_times/t0 +1);
    
    pause_end_indices = round(pause_end_times/t0 +1);
    
    pause_start_indices = pause_start_indices-(T+1);
    
    pause_start_indices(pause_start_indices < 1) = 1;
    
    bool = true([1, len]);
    
    for j=1:length(pause_start_indices)
        
        bool(pause_start_indices(j):pause_end_indices(j)) = false(1);
        
    end
    
    bool = bool(1:length(Dist_array));
    
    Dist_array = Dist_array(bool);
        
    Dists = [Dists, Dist_array];
end

velocities = Dists/(T*t0);

centers = Vmin:binwidth:Vmax;

bootstat = bootstrp(nboot,@ (velocities) V_hist(velocities, centers, binwidth), velocities);

Mean_V_hist=zeros([1, length(centers)]);

Err=zeros([1, length(centers)]);

for i=1:length(centers)
    
    Mean_V_hist(:,i)=mean(bootstat(:,i));
    
    Err(:,i)=std(bootstat(:,i));
end

centers = centers(2:length(centers)-1);

Mean_V_hist = Mean_V_hist(2:length(Mean_V_hist)-1);

Err = Err(2:length(Err)-1);

V_hist_dat = [centers', Mean_V_hist', Err', [binwidth, length(velocities), zeros([1, length(centers)-2])]'];

figure()

ax = axes();

PL1 = errorbar(centers, Mean_V_hist, Err, 'rO', 'MarkerFacecolor', 'r');

set(ax,'XScale','linear','YScale','linear')
    
xlabel('Pause Free Velocity (bp/s)')

ylabel('PDF')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'T = ' num2str(T) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

saveas(PL1, [str_V_pdf, '.jpg'], 'jpg')

saveas(PL1, [str_V_pdf, '.fig'], 'fig')

save([str_V_pdf, '.txt'], 'V_hist_dat', '-ascii')

save([str_V_dat, '.txt'], 'velocities', '-ascii')

