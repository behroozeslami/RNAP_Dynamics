W = 12;

Ts = (2*W+1);

T = Ts;

k = 0;

Nw = 10;

%Pause_threshold =  100.64;
Pause_threshold =  10000;

%Pause_threshold_V = 37.88;

Pause_threshold_V = 100.64;

nboot=100;

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

t0 = 1/cam_freq;

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

mainfoldername = 'Export_AF_7.5pN_1mM_NTP_LP3';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

noise_V_pdf_folder = [output_foldername slash 'Noise_Velocity_pdf' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

if ~(exist(noise_V_pdf_folder,'dir')==7)
    
    mkdir(noise_V_pdf_folder)
    
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

str_noise_V_pdf = [noise_V_pdf_folder slash 'NoiseVelocity_PDF' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise_V_dat = [noise_V_pdf_folder slash 'NoiseVelocity' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise = [noise_V_pdf_folder slash 'NoiseSample' '_' 'Ts='  num2str(2*W+1) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_noise_small = [noise_V_pdf_folder slash 'NoiseSample_small' '_' 'Ts='  num2str(2*W+1) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '.mat'];

x_noise = [];

%sigma_max = 20;

sigma_max = Inf;

 for i=1:length(selected_trace_number)
    
    filename = ['RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = selected_l0_array(i);
    
    factor = l0/l0_WLC;

    x_array = factor*trace_dat(:,8);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);
    
    [Tc_array] = find_crossing_time(t_array, xs_array,Nw);
    
    Tc_array = correct_crossing_time(Tc_array, Pause_threshold, t0);
    
    [DWT_array]=find_Dwell_time(Tc_array);
    
    DWT_indices = 1:length(DWT_array);

    pause_start_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V));
    
    pause_end_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V)+1);
    
    pause_start_indices = round(pause_start_times/t0 +1);
    
    pause_end_indices = round(pause_end_times/t0 +1);
    
    for n=1:length(pause_start_indices)
        
        pause_start = pause_start_indices(n);
        
        pause_end = pause_end_indices(n);
        
        x_chunk = x_array(pause_start:pause_end);
        
        x_chunk = x_chunk(Ts:length(x_chunk)-Ts);
        
        x_chunk = x_chunk - mean(x_chunk);
        
        xs_chunk = xs_array(pause_start:pause_end);
        
        xs_chunk = xs_chunk(Ts:length(xs_chunk)-Ts);
        
        xs_chunk = xs_chunk - mean(xs_chunk);
        
        len = length(xs_chunk);
        
        if std(x_chunk) < sigma_max
            
            x_noise = [x_noise, transpose(x_chunk)];
            
        else 
            
            disp(selected_trace_number(i))
            
        end
    end
        
 end
 
sigma_max = 2*std(x_noise);
 
Dists = [];

x_noise = [];

xs_noise = [];

t_borders = []; 

for i=1:length(selected_trace_number)
    
    filename = ['RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    str_load = [directory slash filename];

    trace_dat = load(str_load , '-ascii');
    
    t_array = trace_dat(:,1);
    
    l0 = selected_l0_array(i);
    
    factor = l0/l0_WLC;

    x_array = factor*trace_dat(:,8);
    
    [xs_array]=Smooth_Trace_SG(x_array,W, k);
    
    [Tc_array] = find_crossing_time(t_array, xs_array,Nw);
    
    Tc_array = correct_crossing_time(Tc_array, Pause_threshold, t0);
    
    [DWT_array]=find_Dwell_time(Tc_array);
    
    DWT_indices = 1:length(DWT_array);

    pause_start_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V));
    
    pause_end_times = Tc_array(DWT_indices(DWT_array - t0/2 > Pause_threshold_V)+1);
    
    pause_start_indices = round(pause_start_times/t0 +1);
    
    pause_end_indices = round(pause_end_times/t0 +1);
    
    for n=1:length(pause_start_indices)
        
        pause_start = pause_start_indices(n);
        
        pause_end = pause_end_indices(n);
        
        x_chunk = x_array(pause_start:pause_end);
        
        x_chunk = x_chunk(Ts:length(x_chunk)-Ts);
        
        x_chunk = x_chunk - mean(x_chunk);
        
        xs_chunk = xs_array(pause_start:pause_end);
        
        xs_chunk = xs_chunk(Ts:length(xs_chunk)-Ts);
        
        xs_chunk = xs_chunk - mean(xs_chunk);
        
        len = length(xs_chunk);
        
        if std(x_chunk) < sigma_max
    
            Dist_array = transpose(xs_chunk(T+1:len) - xs_chunk(1:len-T));
        
            Dists = [Dists, Dist_array];
        
            x_noise = [x_noise, transpose(x_chunk)];
        
            xs_noise = [xs_noise, transpose(xs_chunk)];
        
            t_borders = [t_borders, t0*(length(x_noise)-1)];
            
        end
    end
        
end

t_all = t0*(0:length(x_noise)-1);

Pause_number_array = zeros([1, length(t_all)]);

I1 = 1;

for j=1:length(t_borders)
    
    tb = t_borders(j);
    
    I2 = round(tb/t0)+1;
    
    Pause_number_array(I1:I2) = j;
    
    I1 = I2 + 1;
end

noise_velocities = Dists/(T*t0);

Sym_noise_velocities = [noise_velocities(noise_velocities < 0), -noise_velocities(noise_velocities < 0)];

mu = mean(Sym_noise_velocities);

STD = std(Sym_noise_velocities);

noise_V_dat = Sym_noise_velocities';

%noise_V_dat = noise_velocities';

noise_dat = [x_noise', xs_noise', Pause_number_array'];

noise_dat_small = x_noise';

binwidth = 1;

centers = -60:binwidth:60;

noise_V_hist = V_hist(Sym_noise_velocities, centers, binwidth);

orig_noise_V_hist = V_hist(noise_velocities, centers, binwidth);

centers = centers(2:length(centers)-1);

noise_V_hist = noise_V_hist(2:length(noise_V_hist)-1);

orig_noise_V_hist = orig_noise_V_hist(2:length(orig_noise_V_hist)-1);

gfit = normpdf(centers,mu,STD);

figure()

ax = axes();

PL1 = plot(centers, orig_noise_V_hist, 'bO');

hold on

plot(centers, noise_V_hist, 'rO', 'MarkerFacecolor', 'r');

plot(centers, gfit, 'g', 'LineWidth', 2)

hold off

set(ax,'XScale','linear','YScale','log')
    
xlabel('Noise Velocity (bp/s)')

ylabel('PDF')

legend('Before symmetrization', 'After symmetrization', 'Gaussian fit', 'Location','south')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'T = ' num2str(T) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

figure()

PL2 = plot(t_all, x_noise, 'g')

hold on

plot(t_all, xs_noise, 'r')

hold off

xlabel('Time (s)')

ylabel('Position (bp)')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

saveas(PL1, [str_noise_V_pdf, '.jpg'], 'jpg')

saveas(PL1, [str_noise_V_pdf, '.fig'], 'fig')

save([str_noise_V_dat, '.txt'], 'noise_V_dat', '-ascii')

saveas(PL2, [str_noise, '.jpg'], 'jpg')

%saveas(PL2, [str_noise, '.fig'], 'fig')

save([str_noise, '.txt'], 'noise_dat', '-ascii')

save(str_noise_small, 'noise_dat_small', '-mat')


