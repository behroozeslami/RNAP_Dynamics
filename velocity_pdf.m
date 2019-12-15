W = 12;

Ts = (2*W+1);

T = Ts;

k = 0;

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

mainfoldername = 'EXPORT_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

V_pdf_folder = [output_foldername slash 'Velocity_pdf' '_' 'Ts=' '_' num2str(2*W+1) '_' 'T=' num2str(T)];

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

str_V_pdf = [V_pdf_folder slash 'Velocity_PDF' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

str_negative_V = [V_pdf_folder slash 'negative_Velocities' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)];

str_noise_V_pdf = [V_pdf_folder slash 'Noise_Velocity_PDF(from_necative_V)' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max)]; 

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
    
    xs_array = xs_array((2*W+1):(length(t_array)-(2*W+1)));
    
    t_array = t_array((2*W+1):(length(t_array)-(2*W+1)));
    
    t_array = t_array - t_array(1);

    len = length(xs_array);
    
    Dist_array = transpose(xs_array(T+1:len) - xs_array(1:len-T));
        
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

figure()

ax = axes();

%PL1 = errorbar(centers, Mean_V_hist, Err, 'rO', 'MarkerFacecolor', 'r');

PL1 = plot(centers, Mean_V_hist, 'rO', 'MarkerFacecolor', 'r');

set(ax,'XScale','linear','YScale','log')
    
xlabel('Velocity (bp/s)')

ylabel('PDF')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'T = ' num2str(T) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

saveas(PL1, [str_V_pdf, '.jpg'], 'jpg')

saveas(PL1, [str_V_pdf, '.fig'], 'fig')

negative_velocities = transpose(velocities(velocities < 0));

noise_velocities = [negative_velocities', -negative_velocities'];

mu = mean(noise_velocities);

STD = std(noise_velocities);

binwidth = 1;

centers = -80:binwidth:80;

noise_V_hist = V_hist(noise_velocities, centers, binwidth);

centers = centers(2:length(centers)-1);

noise_V_hist = noise_V_hist(2:length(noise_V_hist)-1);

gfit = normpdf(centers,mu,STD);

figure()

ax = axes();

PL2 = plot(centers, noise_V_hist, 'rO', 'MarkerFacecolor', 'r');

hold on

plot(centers, gfit, 'g', 'LineWidth', 2)

hold off

set(ax,'XScale','linear','YScale','log')
    
xlabel('Velocity (bp/s)')

ylabel('PDF')

str_title = ['Ts = ' ' ' num2str(2*W+1) ' ' 'T = ' num2str(T) ' ' 'l0 = ' num2str(l0_min) '-' num2str(l0_max) ' ' 'Seidel = ' num2str(Seidel_min) '-' num2str(Seidel_max)];

title(str_title)

saveas(PL2, [str_noise_V_pdf, '.jpg'], 'jpg')

saveas(PL2, [str_noise_V_pdf, '.fig'], 'fig')

save([str_negative_V, '.txt'], 'negative_velocities', '-ascii')




