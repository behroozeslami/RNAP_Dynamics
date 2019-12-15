function out = get_DWTpdf_params(DWT_array, Ti_peak, Tf_peak)

peak_DWT = DWT_array((DWT_array > Ti_peak) & (DWT_array < Tf_peak));

Pause_DWT = DWT_array((DWT_array > Tf_peak));

log_peak_DWT = log(peak_DWT);

mu = mean(log_peak_DWT);

sigma_LG = std(log_peak_DWT);

peak_pos = exp(mu - sigma_LG^2);

P = length(Pause_DWT)/length(DWT_array);

out = [peak_pos, sigma_LG, P];