function V_Mat = bootstrap_V(output_foldername, Nw, W, T, t0, Pause_threshold_V,Vmin_fit, Vmax_fit, nboot)

% selection

l0_min = 0.13;

l0_max = 0.46;

Seidel_min = 0.0;

Seidel_max = 0.5;

slash = '/';

str_out = [output_foldername '/V_Bst' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T) '_' 'limits=' num2str(Vmin_fit) '-' num2str(Vmax_fit) '.mat'];

V_pdf_folder = ['../' output_foldername slash 'Velocity_pdf_NoLongPause_v2' '_' 'Ts=' num2str(2*W+1) '_' 'T=' num2str(T) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's'];

str_V_dat = [V_pdf_folder slash 'PauseFreeVelocity' '_' 'Ts='  num2str(2*W+1) '_' 'T=' num2str(T) '_' 'l0_range=' num2str(l0_min) '-' num2str(l0_max) '_' 'Seidel_range=' num2str(Seidel_min) '-' num2str(Seidel_max) '_' 'Nw=' num2str(Nw) '_' 'Pause_threshold_V=' num2str(Pause_threshold_V) 's' '.txt'];

if (exist(str_out, 'file') == 2)
    
    data = load(str_out, '-mat');
    
    V_Mat = data. V_Mat;
    
else
    
    if ~(exist(output_foldername,'dir')==7)
    
        mkdir(output_foldername)
    
    end
    
    Velocities = load(str_V_dat, '-ascii');
    
    Velocities = Velocities(Velocities > Vmin_fit & Velocities < Vmax_fit);
    
    V_Mat = bootstrp(nboot, @(x) x, Velocities);
    
    save(str_out, 'V_Mat', '-mat')
    
end