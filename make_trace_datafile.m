force = 7.5;

cam_freq = 25;
%%%%%% IMPORTANT
%cam_freq = 50;

mainfoldername = 'heterogeneity';

subfolder = 'Fast';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end


if ~(exist(output_foldername,'dir')==7)
    
    mkdir(output_foldername)
    
end

directory = ['..' slash 'Richard_Data' slash mainfoldername slash subfolder];

if strcmp('',subfolder)
    directory = ['..' slash 'Richard_Data' slash mainfoldername];
end

WLC=xlsread(['..' slash 'Richard_Data' slash 'WLC.xlsx']);

force_array = WLC(:,1);

bp_step_array = WLC(:,2);

bp_step_Modelfit = mean(bp_step_array(abs(force_array - force) == min(abs(force_array - force))));

info = dir(directory);

filenames = {info.name};

filenames = filenames(3:length(filenames));

trace_num = length(filenames);

trace_ids = zeros([trace_num, 1]);

Seidel = zeros([trace_num, 1]);

bp_step_Exp = zeros([trace_num, 1]);

bp_step_WLC = bp_step_Modelfit * ones([trace_num, 1]);

for n=1:trace_num
    
    filename = filenames(n);
    
    name = strtok(filename, '.');
    
    id = strtok(name, 'RNAP_');
    
    trace_ids(n,1) = str2num(id{1});
    
    str_load = [directory slash filename{1}];

    dat = load(str_load , '-ascii');
    
    bp_step_Exp(n,1) = dat(1,10);
    
    Seidel(n,1) = dat(2,9);
    
end

output_data = [trace_ids, Seidel, bp_step_Exp, bp_step_WLC];

outfile_str_txt = [output_foldername slash 'Trace_data' '.txt'];

save(outfile_str_txt, 'output_data', '-ascii')

fid = fopen( [output_foldername '/' 'Blacklist.txt'], 'wt' );

fclose(fid);

fid = fopen( [output_foldername '/' 'Blacklist_High_Velocity.txt'], 'wt' );

fclose(fid);

fid = fopen( [output_foldername '/' 'Blacklist_Low_Velocity.txt'], 'wt' );

fclose(fid);