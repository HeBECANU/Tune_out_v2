%% Extracting response phasor - kind of ad hoc at the moment
function signed_signal = phasor_projection(sig,set_freqs,f_thresh)
% INPUTS
%   set_freqs: list of UNIQUE setpoint frequencies
%   sig: Nx2 array of (setpoint, phasor) values
%   f_thresh: Max frequency, below which all phasors are averaged to find
%       phase of modulation response 
% OUTPUTS
%   list of form [f_setpoint, signed_magnitude_response]
    phasor_mask = find(set_freqs < f_thresh);
    phasor_data = sig(phasor_mask,2);
    phasor_mean = mean(phasor_data);
    phasor_std_re = std(real(phasor_data));
    phasor_std_im = std(imag(phasor_data));
    phasor_vector = [real(phasor_mean),imag(phasor_mean)];
    sig_vector = [real(sig(:,2)),imag(sig(:,2))];
    phasor_mat = repmat(phasor_vector, size(sig_vector,1),1);
    phasor_proj = dot(phasor_mat,sig_vector,2);
    phasor_sign = sign(phasor_proj);
    signed_signal = [sig(:,1),abs(sig(:,2)).*phasor_sign];
   
end